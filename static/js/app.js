const { useState, useEffect, useRef } = React;

// Main App Component
function App() {
    const [proteinData, setProteinData] = useState(null);
    const [loading, setLoading] = useState(false);
    const [sequence, setSequence] = useState('');
    const [visualizationStyle, setVisualizationStyle] = useState('sphere');
    const [colorScheme, setColorScheme] = useState('spectrum');
    const [chatMessage, setChatMessage] = useState('');
    const [chatResponse, setChatResponse] = useState('');
    const [simLoading, setSimLoading] = useState(false);
    const [simResults, setSimResults] = useState([]);
    const [simProps, setSimProps] = useState(null);
    
    // Docking-related state
    const [ligandFile, setLigandFile] = useState(null);
    const [ligandData, setLigandData] = useState(null);
    const [dockingResults, setDockingResults] = useState(null);
    const [dockingLoading, setDockingLoading] = useState(false);
    const [showLigandInput, setShowLigandInput] = useState(false);
    const [viewMode, setViewMode] = useState('protein'); // 'protein', 'ligand', 'complex'
    const [originalProteinData, setOriginalProteinData] = useState(null);

    // Load default protein on mount
    useEffect(() => {
        loadDefaultProtein();
    }, []);

    const loadDefaultProtein = async () => {
        setLoading(true);
        try {
            // Load hemoglobin sample
            const response = await fetch('/api/sample/sample_data/hemoglobin.pdb');
            const data = await response.json();

            if (data.success) {
                setProteinData(data.structure);
            } else {
                // Fallback to mock data if sample not available
                const mockProteinData = {
                    header: {
                        title: 'HEMOGLOBIN',
                        resolution: 2.1,
                        experiment_type: 'X-RAY DIFFRACTION',
                        organism: 'Homo sapiens'
                    },
                    atoms: generateMockAtoms(),
                    stats: {
                        atom_count: 4532,
                        chain_count: 4,
                        residue_count: 574,
                        bond_count: 4621
                    }
                };
                setProteinData(mockProteinData);
            }
        } catch (error) {
            console.error('Error loading default protein:', error);
            // Fallback to mock data
            const mockProteinData = {
                header: {
                    title: 'HEMOGLOBIN',
                    resolution: 2.1,
                    experiment_type: 'X-RAY DIFFRACTION',
                    organism: 'Homo sapiens'
                },
                atoms: generateMockAtoms(),
                stats: {
                    atom_count: 4532,
                    chain_count: 4,
                    residue_count: 574,
                    bond_count: 4621
                }
            };
            setProteinData(mockProteinData);
        } finally {
            setLoading(false);
        }
    };

    const generateMockAtoms = () => {
        const atoms = [];
        for (let i = 0; i < 100; i++) {
            atoms.push({
                serial: i + 1,
                name: i % 4 === 0 ? 'CA' : i % 4 === 1 ? 'C' : i % 4 === 2 ? 'N' : 'O',
                res_name: 'ALA',
                chain: 'A',
                res_seq: Math.floor(i / 4) + 1,
                x: (Math.random() - 0.5) * 20,
                y: (Math.random() - 0.5) * 20,
                z: (Math.random() - 0.5) * 20,
                element: i % 4 === 0 ? 'C' : i % 4 === 1 ? 'C' : i % 4 === 2 ? 'N' : 'O',
                occupancy: 1.0,
                temp_factor: 20.0
            });
        }
        return atoms;
    };

    const handleGenerateStructure = async () => {
    if (!sequence.trim()) return;

    setLoading(true);
    try {
        const res = await fetch('/api/sequence', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sequence })
        });
        const data = await res.json();
        if (data.success) {
            setProteinData(data.structure);
        } else {
            setProteinData(null);
            alert(data.error || 'Failed to generate structure.');
        }
    } catch (error) {
        console.error('Error generating structure:', error);
        setProteinData(null);
        alert('Error generating structure.');
    } finally {
        setLoading(false);
    }
};

    const handleFindSimilar = async () => {
        if (!sequence.trim()) return;
        if (sequence.trim().length < 10) {
            alert('Sequence too short for similarity search (min 10 amino acids).');
            return;
        }
        setSimLoading(true);
        try {
            const res = await fetch('/api/similarity', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ sequence })
            });
            const data = await res.json();
            if (data.success) {
                setSimResults(data.hits || []);
                setSimProps(data.top_properties || null);
            } else {
                setSimResults([]);
                setSimProps(null);
                alert(data.error || 'Similarity search failed.');
            }
        } catch (e) {
            console.error('Similarity search error:', e);
            setSimResults([]);
            setSimProps(null);
            alert('Similarity search error.');
        } finally {
            setSimLoading(false);
        }
    };

    const handleChatSend = async () => {
        if (!chatMessage.trim()) return;

        // Send message to backend and get response, include current sequence and simProps as context hints
        const res = await fetch('/api/chat', {
            method: 'POST',
            headers: {'Content-Type': 'application/json' },
            body: JSON.stringify({ 
                message: chatMessage,
                sequence: sequence || '',
                simProps: simProps || null
            })
        });
    const data = await res.json();

    setChatResponse((data.context ? `Context: ${data.context}\n\n` : '') + (data.response || ''));
        setChatMessage('');
    };

    const handleLigandUpload = async (file) => {
        if (!file) return;
        
        const fileExtension = file.name.split('.').pop().toLowerCase();
        if (!['sdf', 'pdb', 'mol', 'mol2'].includes(fileExtension)) {
            alert('Please upload a valid ligand file (SDF, PDB, MOL, MOL2)');
            return;
        }

        setLoading(true);
        try {
            const formData = new FormData();
            formData.append('ligand_file', file);
            
            const res = await fetch('/api/ligand/upload', {
                method: 'POST',
                body: formData
            });
            const data = await res.json();
            
            if (data.success) {
                setLigandData(data.structure);
                setLigandFile(file);
                // Store original protein data if not already stored
                if (!originalProteinData && proteinData) {
                    setOriginalProteinData(proteinData);
                }
                // Automatically switch to ligand view to show the uploaded structure
                setViewMode('ligand');
                const ligandForViewing = {
                    ...data.structure,
                    _isLigand: true
                };
                setProteinData(ligandForViewing);
                alert('Ligand uploaded successfully! Now viewing ligand structure.');
                setShowLigandInput(false);
            } else {
                alert(data.error || 'Failed to upload ligand');
            }
        } catch (error) {
            console.error('Error uploading ligand:', error);
            alert('Error uploading ligand');
        } finally {
            setLoading(false);
        }
    };

    const handleDocking = async () => {
        if (!proteinData) {
            alert('Please ensure protein structure is loaded');
            return;
        }
        if (!ligandFile) {
            alert('Please upload a ligand file first');
            return;
        }

        setDockingLoading(true);
        try {
            const formData = new FormData();
            // Create a temporary protein file from current data
            const proteinBlob = new Blob([convertToPDBFormat(proteinData)], { type: 'text/plain' });
            formData.append('target_file', proteinBlob, 'protein.pdb');
            formData.append('ligand_file', ligandFile);
            formData.append('num_poses', '5');
            formData.append('seed', '42');
            
            const res = await fetch('/api/dock', {
                method: 'POST',
                body: formData
            });
            const data = await res.json();
            
            if (data.success) {
                setDockingResults(data);
                // Store original protein data for view switching
                if (!originalProteinData) {
                    setOriginalProteinData(proteinData);
                }
                // Load the best pose into the viewer
                if (data.best_pose_structure) {
                    setProteinData(data.best_pose_structure);
                    setViewMode('complex');
                }
                alert(`Docking completed! ${data.results.length} poses generated. Best pose loaded.`);
            } else {
                alert(data.error || 'Docking failed');
            }
        } catch (error) {
            console.error('Error performing docking:', error);
            alert('Error performing docking');
        } finally {
            setDockingLoading(false);
        }
    };

    const handleViewModeChange = (mode) => {
        setViewMode(mode);
        
        switch (mode) {
            case 'protein':
                if (originalProteinData) {
                    setProteinData(originalProteinData);
                }
                break;
            case 'ligand':
                if (ligandData) {
                    // Create a copy of ligand data optimized for small molecule viewing
                    const ligandForViewing = {
                        ...ligandData,
                        _isLigand: true  // Flag to help the viewer handle this as a small molecule
                    };
                    setProteinData(ligandForViewing);
                }
                break;
            case 'complex':
                if (dockingResults && dockingResults.best_pose_structure) {
                    setProteinData(dockingResults.best_pose_structure);
                }
                break;
        }
    };

    // Helper function to convert structure data to PDB format
    const convertToPDBFormat = (data) => {
        let pdbLines = ['TITLE     ' + (data.header?.title || 'PROTEIN STRUCTURE')];
        
        data.atoms.forEach(atom => {
            // Use HETATM for ligands or if the original atom was a HETATM
            const recordType = (data._isLigand || atom.is_hetatm) ? 'HETATM' : 'ATOM  ';
            
            // Ensure numeric fields are properly handled
            const serial = parseInt(atom.serial) || 1;
            const resSeq = parseInt(atom.res_seq) || 1;
            const x = parseFloat(atom.x) || 0.0;
            const y = parseFloat(atom.y) || 0.0;
            const z = parseFloat(atom.z) || 0.0;
            const occupancy = parseFloat(atom.occupancy) || 1.0;
            const tempFactor = parseFloat(atom.temp_factor) || 20.0;
            
            const line = [
                recordType,
                serial.toString().padStart(5),
                '  ',
                atom.name.padEnd(4),
                ' ',
                atom.res_name.padEnd(3),
                ' ',
                atom.chain || 'A',
                resSeq.toString().padStart(4),
                '    ',
                x.toFixed(3).padStart(8),
                y.toFixed(3).padStart(8),
                z.toFixed(3).padStart(8),
                occupancy.toFixed(2).padStart(6),
                tempFactor.toFixed(2).padStart(6),
                '          ',
                atom.element.padStart(2)
            ].join('');
            pdbLines.push(line);
        });
        
        pdbLines.push('END');
        return pdbLines.join('\n');
    };

    const handleKeyPress = (e, handler) => {
        if (e.key === 'Enter') {
            handler();
        }
    };

    return (
        <div>
            <TopBar 
                sequence={sequence}
                setSequence={setSequence}
                onGenerate={handleGenerateStructure}
                onFindSimilar={handleFindSimilar}
                loading={loading}
                simLoading={simLoading}
                // Ligand and docking props
                ligandFile={ligandFile}
                onLigandUpload={handleLigandUpload}
                onDocking={handleDocking}
                dockingLoading={dockingLoading}
                showLigandInput={showLigandInput}
                setShowLigandInput={setShowLigandInput}
                ligandData={ligandData}
                // View mode props
                viewMode={viewMode}
                onViewModeChange={handleViewModeChange}
                dockingResults={dockingResults}
            />
            
            <div className="main-container">
                <div className="viewer-section">
                    <ProteinViewer
                        proteinData={proteinData}
                        visualizationStyle={visualizationStyle}
                        colorScheme={colorScheme}
                        loading={loading}
                    />
                </div>
                
                <Sidebar
                    proteinData={proteinData}
                    visualizationStyle={visualizationStyle}
                    setVisualizationStyle={setVisualizationStyle}
                    colorScheme={colorScheme}
                    setColorScheme={setColorScheme}
                    simResults={simResults}
                    simProps={simProps}
                    simLoading={simLoading}
                />
            </div>
            
            <ChatSection
                message={chatMessage}
                setMessage={setChatMessage}
                onSend={handleChatSend}
                response={chatResponse}
                onKeyPress={handleKeyPress}
            />
        </div>
    );
}

// Top Bar Component
function TopBar({ 
    sequence, setSequence, onGenerate, onFindSimilar, loading, simLoading,
    ligandFile, onLigandUpload, onDocking, dockingLoading,
    showLigandInput, setShowLigandInput, ligandData,
    viewMode, onViewModeChange, dockingResults 
}) {
    const ligandInputRef = React.useRef();
    return (
        <div className="top-bar">
            <div className="top-bar-main">
                <div className="app-title">
                    <i className="fas fa-dna"></i>
                    Protein Structure Viewer
                </div>
                
                <div className="sequence-input-container">
                    <input
                        type="text"
                        className="sequence-input"
                        placeholder="Enter amino acid sequence (e.g., MVLSPADKTNVKAAW...)"
                        value={sequence}
                        onChange={(e) => setSequence(e.target.value)}
                        disabled={loading}
                    />
                    
                    <button
                        className="add-ligand-btn"
                        onClick={() => setShowLigandInput(!showLigandInput)}
                        title="Add ligand for docking"
                        style={{ marginLeft: '0.5rem' }}
                    >
                        <i className="fas fa-plus"></i>
                    </button>
                </div>
                
                <div className="button-container">
                    <button 
                        className="generate-btn" 
                        onClick={onGenerate}
                        disabled={loading}
                    >
                        {loading ? (
                            <>
                                <div className="loading-spinner" style={{ width: '1rem', height: '1rem', margin: 0 }}></div>
                                Generating...
                            </>
                        ) : (
                            <>
                                <i className="fas fa-play"></i>
                                Generate
                            </>
                        )}
                    </button>

                    <button
                        className="generate-btn"
                        onClick={onFindSimilar}
                        disabled={simLoading || !sequence.trim()}
                        style={{ marginLeft: '0.5rem' }}
                    >
                        {simLoading ? (
                            <>
                                <div className="loading-spinner" style={{ width: '1rem', height: '1rem', margin: 0 }}></div>
                                Searching...
                            </>
                        ) : (
                            <>
                                <i className="fas fa-search"></i>
                                Find Similar
                            </>
                        )}
                    </button>

                    {ligandFile && (
                        <button
                            className="generate-btn dock-btn"
                            onClick={onDocking}
                            disabled={dockingLoading || !ligandFile}
                            style={{ marginLeft: '0.5rem', backgroundColor: '#059669' }}
                        >
                            {dockingLoading ? (
                                <>
                                    <div className="loading-spinner" style={{ width: '1rem', height: '1rem', margin: 0 }}></div>
                                    Docking...
                                </>
                            ) : (
                                <>
                                    <i className="fas fa-link"></i>
                                    Dock Proteins
                                </>
                            )}
                        </button>
                    )}
                </div>
            </div>
            
            {showLigandInput && (
                <div className="ligand-input-container">
                    <input
                        ref={ligandInputRef}
                        type="file"
                        accept=".sdf,.pdb,.mol,.mol2"
                        onChange={(e) => onLigandUpload(e.target.files[0])}
                        disabled={loading}
                        style={{ display: 'none' }}
                    />
                    <div className="file-upload-wrapper">
                        <button 
                            className="file-upload-btn" 
                            onClick={() => ligandInputRef.current?.click()}
                            disabled={loading}
                        >
                            <i className="fas fa-upload"></i>
                            {ligandFile ? ligandFile.name : 'Choose Ligand File (SDF, PDB, MOL, MOL2)'}
                        </button>
                        {ligandFile && (
                            <span className="file-status">
                                <i className="fas fa-check-circle" style={{ color: '#10b981', marginLeft: '0.5rem' }}></i>
                                Uploaded
                            </span>
                        )}
                    </div>
                </div>
            )}

            {(ligandData || dockingResults) && (
                <div className="view-mode-container">
                    <span className="view-mode-label">View:</span>
                    <div className="view-mode-buttons">
                        <button
                            className={`view-mode-btn ${viewMode === 'protein' ? 'active' : ''}`}
                            onClick={() => onViewModeChange('protein')}
                        >
                            <i className="fas fa-dna"></i>
                            Protein
                        </button>
                        <button
                            className={`view-mode-btn ${viewMode === 'ligand' ? 'active' : ''}`}
                            onClick={() => onViewModeChange('ligand')}
                            disabled={!ligandData}
                        >
                            <i className="fas fa-circle"></i>
                            Ligand
                        </button>
                        {dockingResults && (
                            <button
                                className={`view-mode-btn ${viewMode === 'complex' ? 'active' : ''}`}
                                onClick={() => onViewModeChange('complex')}
                            >
                                <i className="fas fa-link"></i>
                                Complex
                            </button>
                        )}
                    </div>
                </div>
            )}
        </div>
    );
}

// Protein Viewer Component
function ProteinViewer({ proteinData, visualizationStyle, colorScheme, loading }) {
    const viewerRef = useRef();
    const molViewerRef = useRef();

    useEffect(() => {
        if (viewerRef.current && !molViewerRef.current && window.$3Dmol) {
            molViewerRef.current = window.$3Dmol.createViewer(viewerRef.current, {
                defaultcolors: window.$3Dmol.rasmolElementColors
            });
            molViewerRef.current.setBackgroundColor('#1e1b4b');
        }
    }, []);

    useEffect(() => {
        if (molViewerRef.current && proteinData) {
            loadProteinData();
        }
    }, [proteinData]);

    useEffect(() => {
        if (molViewerRef.current && proteinData) {
            applyVisualizationStyle();
        }
    }, [visualizationStyle, colorScheme]);

    const loadProteinData = () => {
        if (!molViewerRef.current || !proteinData) return;

        molViewerRef.current.clear();

        try {
            const pdbData = convertToPDBFormat(proteinData);
            const model = molViewerRef.current.addModel(pdbData, 'pdb');
            applyVisualizationStyle();
            molViewerRef.current.zoomTo();
            molViewerRef.current.render();
        } catch (error) {
            console.error('Error loading protein:', error);
        }
    };

    const convertToPDBFormat = (data) => {
        let pdbLines = ['TITLE     ' + (data.header?.title || 'PROTEIN STRUCTURE')];
        
        data.atoms.forEach(atom => {
            // Use HETATM for ligands or if the original atom was a HETATM
            const recordType = (data._isLigand || atom.is_hetatm) ? 'HETATM' : 'ATOM  ';
            
            // Ensure numeric fields are properly handled
            const serial = parseInt(atom.serial) || 1;
            const resSeq = parseInt(atom.res_seq) || 1;
            const x = parseFloat(atom.x) || 0.0;
            const y = parseFloat(atom.y) || 0.0;
            const z = parseFloat(atom.z) || 0.0;
            const occupancy = parseFloat(atom.occupancy) || 1.0;
            const tempFactor = parseFloat(atom.temp_factor) || 20.0;
            
            const line = [
                recordType,
                serial.toString().padStart(5),
                '  ',
                atom.name.padEnd(4),
                ' ',
                atom.res_name.padEnd(3),
                ' ',
                atom.chain || 'A',
                resSeq.toString().padStart(4),
                '    ',
                x.toFixed(3).padStart(8),
                y.toFixed(3).padStart(8),
                z.toFixed(3).padStart(8),
                occupancy.toFixed(2).padStart(6),
                tempFactor.toFixed(2).padStart(6),
                '          ',
                atom.element.padStart(2)
            ].join('');
            pdbLines.push(line);
        });
        
        pdbLines.push('END');
        return pdbLines.join('\n');
    };

    const applyVisualizationStyle = () => {
        if (!molViewerRef.current || !proteinData) return;

        molViewerRef.current.setStyle({}, {});

        const scheme = colorScheme === 'element' ? 'default' : colorScheme;
        const isLigand = proteinData._isLigand || (proteinData.atoms && proteinData.atoms.length < 100);

        if (isLigand) {
            // For ligands/small molecules, use ball-and-stick representation
            molViewerRef.current.setStyle({}, {
                stick: {
                    colorscheme: 'default',
                    radius: 0.3
                },
                sphere: {
                    colorscheme: 'default',
                    scale: 0.3
                }
            });
        } else {
            // For proteins, use the selected visualization style
            switch (visualizationStyle) {
                case 'cartoon':
                    molViewerRef.current.setStyle({}, {
                        cartoon: {
                            color: scheme === 'spectrum' ? 'spectrum' : scheme,
                            thickness: 0.4,
                            opacity: 0.9
                        }
                    });
                    break;
                case 'stick':
                    molViewerRef.current.setStyle({}, {
                        stick: {
                            colorscheme: scheme,
                            radius: 0.25
                        }
                    });
                    break;
                case 'sphere':
                    molViewerRef.current.setStyle({}, {
                        sphere: {
                            colorscheme: scheme,
                            scale: 0.8,
                            opacity: 0.9
                        }
                    });
                    break;
                case 'line':
                    molViewerRef.current.setStyle({ atom: 'CA' }, {
                        line: {
                            colorscheme: scheme,
                            linewidth: 2
                        }
                    });
                    break;
            }
        }

        molViewerRef.current.render();
    };

    const handleReset = () => {
        if (molViewerRef.current) {
            molViewerRef.current.zoomTo();
            molViewerRef.current.render();
        }
    };

    const handleFullscreen = () => {
        if (viewerRef.current.requestFullscreen) {
            viewerRef.current.requestFullscreen();
        }
    };

    if (loading) {
        return (
            <div className="protein-viewer">
                <div className="loading">
                    <div className="loading-spinner"></div>
                    <p>Loading protein structure...</p>
                </div>
            </div>
        );
    }

    return (
        <div className="protein-viewer">
            <div 
                ref={viewerRef} 
                style={{ width: '100%', height: '100%' }}
            />
            
            {proteinData && (
                <>
                    <div className="viewer-overlay">
                        <h3>{proteinData.header.title}</h3>
                        <p>{proteinData.stats.atom_count} atoms, {proteinData.stats.chain_count} chains</p>
                    </div>
                    
                    <div className="viewer-controls">
                        <button className="control-btn" onClick={handleReset} title="Reset View">
                            <i className="fas fa-home"></i>
                        </button>
                        <button className="control-btn" onClick={handleFullscreen} title="Fullscreen">
                            <i className="fas fa-expand"></i>
                        </button>
                    </div>
                </>
            )}
            
            {!proteinData && !loading && (
                <div className="loading">
                    <i className="fas fa-dna" style={{ fontSize: '3rem', marginBottom: '1rem', color: '#60a5fa' }}></i>
                    <p>Enter a sequence to generate structure</p>
                </div>
            )}
        </div>
    );
}

// Sidebar Component
function Sidebar({ proteinData, visualizationStyle, setVisualizationStyle, colorScheme, setColorScheme, simResults, simProps, simLoading }) {
    const styles = [
    { value: 'stick', label: 'Stick', icon: 'fas fa-grip-lines' },
    { value: 'sphere', label: 'Spacefill', icon: 'fas fa-dot-circle' },
    { value: 'line', label: 'Backbone', icon: 'fas fa-project-diagram' }
    ];

    // const colors = [
    //     { value: 'spectrum', label: 'Spectrum' },
    //     { value: 'chain', label: 'By Chain' },
    //     { value: 'residue', label: 'By Residue' },
    //     { value: 'element', label: 'By Element' }
    // ];

    return (
        <div className="sidebar">
            <div className="sidebar-section">
                <h3>
                    Visualization Style
                </h3>
                
                <div className="control-group">
                    <label className="control-label">Representation</label>
                    <div className="control-buttons">
                        {styles.map((style) => (
                            <button
                                key={style.value}
                                className={`style-btn ${visualizationStyle === style.value ? 'active' : ''}`}
                                onClick={() => setVisualizationStyle(style.value)}
                            >
                                <i className={style.icon}></i>
                                <span>{style.label}</span>
                            </button>
                        ))}
                    </div>
                </div>
                
                {/* <div className="control-group">
                    <label className="control-label">Color Scheme</label>
                    <select 
                        className="color-select"
                        value={colorScheme} 
                        onChange={(e) => setColorScheme(e.target.value)}
                    >
                        {colors.map((color) => (
                            <option key={color.value} value={color.value}>
                                {color.label}
                            </option>
                        ))}
                    </select>
                </div> */}
            </div>

            {proteinData && (
                <div className="sidebar-section">
                    <h3>
                        <i className="fas fa-info-circle"></i>
                        Protein Information
                    </h3>
                    
                    <div className="info-grid">
                        <div className="info-item">
                            <span className="info-label">Atoms</span>
                            <span className="info-value">{proteinData.stats.atom_count.toLocaleString()}</span>
                        </div>
                        
                        <div className="info-item">
                            <span className="info-label">Chains</span>
                            <span className="info-value">{proteinData.stats.chain_count}</span>
                        </div>
                        
                        <div className="info-item">
                            <span className="info-label">Residues</span>
                            <span className="info-value">{proteinData.stats.residue_count}</span>
                        </div>
                        
                        {proteinData.header.resolution && (
                            <div className="info-item">
                                <span className="info-label">Resolution</span>
                                <span className="info-value">{proteinData.header.resolution} Å</span>
                            </div>
                        )}
                        
                        {proteinData.header.experiment_type && (
                            <div className="info-item">
                                <span className="info-label">Method</span>
                                <span className="info-value">{proteinData.header.experiment_type}</span>
                            </div>
                        )}
                        
                        {proteinData.header.organism && (
                            <div className="info-item">
                                <span className="info-label">Organism</span>
                                <span className="info-value">{proteinData.header.organism}</span>
                            </div>
                        )}
                    </div>
                </div>
            )}

            <div className="sidebar-section">
                <h3>
                    <i className="fas fa-search"></i>
                    Similarity Search
                </h3>
                {simLoading && (
                    <div className="loading" style={{ padding: '0.5rem 0' }}>
                        <div className="loading-spinner"></div>
                        <p>Searching...</p>
                    </div>
                )}
                {!simLoading && (!simResults || simResults.length === 0) && (
                    <p style={{ opacity: 0.7 }}>Run "Find Similar" to see top matches.</p>
                )}

                {simProps && (
                    <div style={{ marginTop: '0.5rem' }}>
                        <h4 style={{ margin: '0.25rem 0' }}>Top-hit properties</h4>
                        <div className="info-grid">
                            <div className="info-item"><span className="info-label">Length</span><span className="info-value">{simProps.length}</span></div>
                            <div className="info-item"><span className="info-label">Mass</span><span className="info-value">{simProps.mass_Da} Da</span></div>
                            <div className="info-item"><span className="info-label">Hydropathy</span><span className="info-value">{simProps.hydropathy_KD}</span></div>
                            <div className="info-item"><span className="info-label">Net charge ~pH7</span><span className="info-value">{simProps.net_charge_pH7_approx}</span></div>
                            <div className="info-item"><span className="info-label">Aromatic</span><span className="info-value">{simProps.aromatic_count}</span></div>
                        </div>
                    </div>
                )}
            </div>
        </div>
    );
}

// Chat Section Component
function ChatSection({ message, setMessage, onSend, response, onKeyPress }) {
    return (
        <div className="chat-section">
            {response && (
                <div className="chat-response" style={{ whiteSpace: 'pre-wrap' }}>
                    <strong>AI:</strong> {response}
                </div>
            )}
            
            <div className="chat-input-container">
                <input
                    type="text"
                    className="chat-input"
                    placeholder="Ask me anything about this protein..."
                    value={message}
                    onChange={(e) => setMessage(e.target.value)}
                    onKeyPress={(e) => onKeyPress(e, onSend)}
                />
                
                <button 
                    className="chat-send-btn" 
                    onClick={onSend}
                    disabled={!message.trim()}
                >
                    <i className="fas fa-paper-plane"></i>
                </button>
            </div>
        </div>
    );
}

// Render the app
ReactDOM.render(<App />, document.getElementById('root'));