"""
Protein similarity search utilities.

Features:
- Load a small database of short protein sequences (FASTA format).
- Local alignment scoring (Smith-Waterman style) with BLOSUM62 and linear gap penalty.
- Top-K search returning percent identity and score.
- Basic property computation for the best hit (length, mass, hydrophobicity, charge).

Notes:
- Designed for short proteins/peptides (< 100 aa). For longer sequences, consider
  an optimized library (e.g., parasail) or BLAST.
- Database path can be configured via environment variable PROTEIN_DB_FASTA.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# 20 standard amino acids
AA = "ARNDCQEGHILKMFPSTWYV"


# BLOSUM62 substitution matrix (integer scores)
# Source: Henikoff & Henikoff (1992); public domain data table.
_BLOSUM62_ROWS = AA
_BLOSUM62 = [
    #   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],  # A
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],  # R
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],  # N
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],  # D
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],  # C
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],  # Q
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],  # E
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],  # G
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],  # H
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],  # I
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],  # L
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],  # K
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],  # M
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],  # F
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],  # P
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],  # S
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],  # T
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],  # W
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],  # Y
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],  # V
]

_B62: Dict[Tuple[str, str], int] = {}
for i, a in enumerate(_BLOSUM62_ROWS):
    for j, b in enumerate(_BLOSUM62_ROWS):
        _B62[(a, b)] = _BLOSUM62[i][j]


# Kyte-Doolittle hydropathy index
KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2
}

# Average residue masses (Daltons) of free amino acids (approximate)
AA_MASS = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}


@dataclass
class Protein:
    id: str
    name: str
    sequence: str
    organism: Optional[str] = None
    function: Optional[str] = None


DB: List[Protein] = []


def is_valid_protein_sequence(seq: str) -> bool:
    if not seq:
        return False
    s = seq.strip().upper()
    return all(c in AA for c in s)


def _score(a: str, b: str) -> int:
    return _B62.get((a, b), -1)


def smith_waterman(a: str, b: str, gap: int = -5) -> Tuple[int, str, str, int]:
    """
    Local alignment with linear gap penalty.
    Returns: (best_score, aligned_a, aligned_b, matches)
    """
    a = a.upper()
    b = b.upper()
    n, m = len(a), len(b)
    # H: score matrix; P: pointer (0 stop, 1 diag, 2 up, 3 left)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[0] * (m + 1) for _ in range(n + 1)]

    best = 0
    bi = bj = 0

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            bjc = b[j - 1]
            s = _score(ai, bjc)
            d = H[i - 1][j - 1] + s
            u = H[i - 1][j] + gap
            l = H[i][j - 1] + gap
            val = max(0, d, u, l)
            H[i][j] = val
            if val == 0:
                P[i][j] = 0
            elif val == d:
                P[i][j] = 1
            elif val == u:
                P[i][j] = 2
            else:
                P[i][j] = 3
            if val > best:
                best = val
                bi, bj = i, j

    # Traceback from (bi,bj)
    aligned_a = []
    aligned_b = []
    matches = 0
    i, j = bi, bj
    while i > 0 and j > 0 and P[i][j] != 0:
        if P[i][j] == 1:
            ai = a[i - 1]
            bjc = b[j - 1]
            aligned_a.append(ai)
            aligned_b.append(bjc)
            if ai == bjc:
                matches += 1
            i -= 1
            j -= 1
        elif P[i][j] == 2:
            aligned_a.append(a[i - 1])
            aligned_b.append('-')
            i -= 1
        else:
            aligned_a.append('-')
            aligned_b.append(b[j - 1])
            j -= 1

    return best, ''.join(reversed(aligned_a)), ''.join(reversed(aligned_b)), matches


def percent_identity(aln_a: str, aln_b: str) -> float:
    length = sum(1 for x, y in zip(aln_a, aln_b) if x != '-' or y != '-')
    if length == 0:
        return 0.0
    matches = sum(1 for x, y in zip(aln_a, aln_b) if x == y and x != '-')
    return 100.0 * matches / length


def load_fasta(path: str) -> List[Protein]:
    proteins: List[Protein] = []
    if not os.path.exists(path):
        return proteins
    with open(path, 'r') as f:
        header = None
        seq_lines: List[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    proteins.append(_protein_from_header_seq(header, ''.join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            proteins.append(_protein_from_header_seq(header, ''.join(seq_lines)))
    return proteins


def _protein_from_header_seq(header: str, seq: str) -> Protein:
    # Header format: id|name|organism|function  (fields after id optional and can be '-')
    parts = header.split('|')
    pid = parts[0].strip()
    name = parts[1].strip() if len(parts) > 1 else pid
    org = parts[2].strip() if len(parts) > 2 and parts[2].strip() != '-' else None
    func = parts[3].strip() if len(parts) > 3 and parts[3].strip() != '-' else None
    return Protein(id=pid, name=name or pid, sequence=seq.strip().upper(), organism=org, function=func)


def init_database(fasta_path: Optional[str] = None) -> int:
    """Load the protein database; returns count loaded."""
    global DB
    if fasta_path is None:
        fasta_path = os.getenv('PROTEIN_DB_FASTA', os.path.join(os.path.dirname(__file__), 'data', 'proteins_demo.fasta'))
    DB = load_fasta(fasta_path)
    # Filter to short proteins
    DB = [p for p in DB if 1 <= len(p.sequence) <= 100 and is_valid_protein_sequence(p.sequence)]
    return len(DB)


@dataclass
class Hit:
    protein: Protein
    score: int
    identity_percent: float
    alignment_length: int
    coverage: float
    aligned_query: str
    aligned_subject: str


def search(query_seq: str, top_k: int = 3, min_coverage: float = 0.3, min_aln_len: int = 12) -> List[Hit]:
    if not DB:
        init_database()
    if not is_valid_protein_sequence(query_seq):
        raise ValueError('Invalid protein sequence: only 20 standard amino acids allowed (ACDEFGHIKLMNPQRSTVWY).')
    q = query_seq.strip().upper()
    hits: List[Hit] = []
    for p in DB:
        score, qa, sb, _ = smith_waterman(q, p.sequence)
        pid = percent_identity(qa, sb)
        aln_len = sum(1 for x, y in zip(qa, sb) if x != '-' or y != '-')
    cov = aln_len / max(1, len(q))
    if aln_len >= min_aln_len and cov >= min_coverage:
        hits.append(Hit(protein=p, score=score, identity_percent=pid, alignment_length=aln_len, coverage=round(cov, 3), aligned_query=qa, aligned_subject=sb))
    # Sort primarily by score, then by coverage, then identity
    hits.sort(key=lambda h: (h.score, h.coverage, h.identity_percent), reverse=True)
    return hits[: top_k]


def compute_properties(seq: str) -> Dict[str, object]:
    s = seq.strip().upper()
    length = len(s)
    mass = sum(AA_MASS.get(c, 0.0) for c in s)
    hydropathy = sum(KD.get(c, 0.0) for c in s) / length if length else 0.0
    comp: Dict[str, int] = {a: 0 for a in AA}
    for c in s:
        if c in comp:
            comp[c] += 1
    # Very rough net charge at pH 7 estimate (H counted half positive)
    pos = s.count('K') + s.count('R') + 0.1 * s.count('H')
    neg = s.count('D') + s.count('E')
    net_charge = pos - neg
    aromatic = s.count('F') + s.count('W') + s.count('Y')
    return {
        'length': length,
        'mass_Da': round(mass, 2),
        'hydropathy_KD': round(hydropathy, 3),
        'net_charge_pH7_approx': round(net_charge, 2),
        'aromatic_count': aromatic,
        'composition': comp,
    }


def summarize_hits(hits: List[Hit]) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for h in hits:
        out.append({
            'id': h.protein.id,
            'name': h.protein.name,
            'organism': h.protein.organism,
            'function': h.protein.function,
            'length': len(h.protein.sequence),
            'score': h.score,
            'identity_percent': round(h.identity_percent, 2),
            'alignment_length': h.alignment_length,
            'coverage': h.coverage,
        })
    return out


if __name__ == '__main__':
    # Simple CLI test
    cnt = init_database()
    print(f"Loaded {cnt} proteins in DB")
    q = os.environ.get('TEST_SEQ', 'MALWMRLLPLLALLALWGPDPAAA')  # insulin signal peptide fragment
    hits = search(q)
    print(summarize_hits(hits))
    top = hits[0] if hits else None
    if top:
        print('Top properties:', compute_properties(top.protein.sequence))
