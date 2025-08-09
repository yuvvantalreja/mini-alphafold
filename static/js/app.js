const { useState, useEffect, useRef } = React;

// Main App Component
function App() {
    const [proteinData, setProteinData] = useState(null);
    const [loading, setLoading] = useState(false);
    const [sequence, setSequence] = useState('');
    const [visualizationStyle, setVisualizationStyle] = useState('cartoon');
    const [colorScheme, setColorScheme] = useState('spectrum');
    const [chatMessage, setChatMessage] = useState('');
    const [chatResponse, setChatResponse] = useState('');

    // Load default protein on mount
    useEffect(() => {
        loadDefaultProtein();
    }, []);

    const loadDefaultProtein = async () => {
        setLoading(true);
        try {
            // Load hemoglobin sample
            const response = await fetch('/api/sample/hemoglobin');
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
            // TODO: Connect to your backend API that generates PDB from sequence
            // For now, simulate the call and reload default protein
            await new Promise(resolve => setTimeout(resolve, 2000));
            await loadDefaultProtein();
        } catch (error) {
            console.error('Error generating structure:', error);
        } finally {
            setLoading(false);
        }
    };

    const handleChatSend = async () => {
        if (!chatMessage.trim()) return;
        
        try {
            // TODO: Connect to your LLM backend
            // For now, always return "Hello" as requested
            setChatResponse('Hello');
            setChatMessage('');
        } catch (error) {
            console.error('Error sending chat message:', error);
        }
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
                loading={loading}
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
function TopBar({ sequence, setSequence, onGenerate, loading }) {
    return (
        <div className="top-bar">
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
            </div>
            
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
        let pdbLines = ['TITLE     ' + (data.header.title || 'PROTEIN STRUCTURE')];
        
        data.atoms.forEach(atom => {
            const line = [
                'ATOM  ',
                atom.serial.toString().padStart(5),
                '  ',
                atom.name.padEnd(4),
                ' ',
                atom.res_name.padEnd(3),
                ' ',
                atom.chain || 'A',
                atom.res_seq.toString().padStart(4),
                '    ',
                atom.x.toFixed(3).padStart(8),
                atom.y.toFixed(3).padStart(8),
                atom.z.toFixed(3).padStart(8),
                atom.occupancy.toFixed(2).padStart(6),
                atom.temp_factor.toFixed(2).padStart(6),
                '          ',
                atom.element.padStart(2)
            ].join('');
            pdbLines.push(line);
        });
        
        pdbLines.push('END');
        return pdbLines.join('\n');
    };

    const applyVisualizationStyle = () => {
        if (!molViewerRef.current) return;

        molViewerRef.current.setStyle({}, {});

        const scheme = colorScheme === 'element' ? 'default' : colorScheme;

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
function Sidebar({ proteinData, visualizationStyle, setVisualizationStyle, colorScheme, setColorScheme }) {
    const styles = [
        { value: 'cartoon', label: 'Cartoon', icon: 'fas fa-ribbon' },
        { value: 'stick', label: 'Stick', icon: 'fas fa-grip-lines' },
        { value: 'sphere', label: 'Spacefill', icon: 'fas fa-dot-circle' },
        { value: 'line', label: 'Backbone', icon: 'fas fa-project-diagram' }
    ];

    const colors = [
        { value: 'spectrum', label: 'Spectrum' },
        { value: 'chain', label: 'By Chain' },
        { value: 'residue', label: 'By Residue' },
        { value: 'element', label: 'By Element' }
    ];

    return (
        <div className="sidebar">
            <div className="sidebar-section">
                <h3>
                    <i className="fas fa-palette"></i>
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
                
                <div className="control-group">
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
                </div>
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
        </div>
    );
}

// Chat Section Component
function ChatSection({ message, setMessage, onSend, response, onKeyPress }) {
    return (
        <div className="chat-section">
            {response && (
                <div className="chat-response">
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