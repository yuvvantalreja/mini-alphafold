const { useState, useEffect, useRef, useCallback } = React;

// Main App Component
function App() {
    const [proteinData, setProteinData] = useState(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    const [success, setSuccess] = useState(null);
    const [visualizationStyle, setVisualizationStyle] = useState('cartoon');
    const [colorScheme, setColorScheme] = useState('spectrum');
    const [samples, setSamples] = useState([]);
    const [uploadedFileName, setUploadedFileName] = useState(null);

    // Load sample proteins on component mount
    useEffect(() => {
        fetchSamples();
    }, []);

    const fetchSamples = async () => {
        try {
            const response = await fetch('/api/samples');
            const data = await response.json();
            setSamples(data.samples || []);
        } catch (err) {
            console.error('Failed to load samples:', err);
        }
    };

    const clearMessages = () => {
        setError(null);
        setSuccess(null);
    };

    const handleFileUpload = async (file) => {
        if (!file) return;

        setLoading(true);
        clearMessages();

        const formData = new FormData();
        formData.append('file', file);

        try {
            const response = await fetch('/api/upload', {
                method: 'POST',
                body: formData
            });

            const data = await response.json();

            if (data.success) {
                setProteinData(data.structure);
                setUploadedFileName(data.filename);
                setSuccess(`Successfully loaded ${data.filename}`);
            } else {
                setError(data.error || 'Failed to upload file');
            }
        } catch (err) {
            setError('Network error: ' + err.message);
        } finally {
            setLoading(false);
        }
    };

    const handleSampleLoad = async (sampleName) => {
        setLoading(true);
        clearMessages();

        try {
            const response = await fetch(`/api/sample/${sampleName}`);
            const data = await response.json();

            if (data.success) {
                setProteinData(data.structure);
                setUploadedFileName(data.filename);
                setSuccess(`Successfully loaded ${data.filename}`);
            } else {
                setError(data.error || 'Failed to load sample');
            }
        } catch (err) {
            setError('Network error: ' + err.message);
        } finally {
            setLoading(false);
        }
    };

    return (
        <div className="app-container">
            <Header 
                proteinData={proteinData}
                uploadedFileName={uploadedFileName}
            />
            
            <Sidebar
                onFileUpload={handleFileUpload}
                proteinData={proteinData}
                visualizationStyle={visualizationStyle}
                setVisualizationStyle={setVisualizationStyle}
                colorScheme={colorScheme}
                setColorScheme={setColorScheme}
                samples={samples}
                onSampleLoad={handleSampleLoad}
                loading={loading}
                error={error}
                success={success}
            />
            
            <ProteinViewer
                proteinData={proteinData}
                visualizationStyle={visualizationStyle}
                colorScheme={colorScheme}
                loading={loading}
            />
        </div>
    );
}

// Header Component
function Header({ proteinData, uploadedFileName }) {
    return (
        <header className="app-header">
            <div className="app-title">
                <i className="fas fa-dna"></i>
                Protein Structure Viewer
            </div>
            <div className="header-controls">
                {uploadedFileName && (
                    <span className="text-sm text-gray-600">
                        <i className="fas fa-file"></i>
                        {uploadedFileName}
                    </span>
                )}
                {proteinData && (
                    <span className="text-sm text-gray-600">
                        {proteinData.stats.atom_count} atoms, {proteinData.stats.chain_count} chains
                    </span>
                )}
            </div>
        </header>
    );
}

// Sidebar Component
function Sidebar({ 
    onFileUpload, 
    proteinData, 
    visualizationStyle, 
    setVisualizationStyle,
    colorScheme,
    setColorScheme,
    samples,
    onSampleLoad,
    loading,
    error,
    success
}) {
    return (
        <div className="sidebar">
            {error && (
                <div className="error">
                    <i className="fas fa-exclamation-triangle"></i>
                    {error}
                </div>
            )}
            
            {success && (
                <div className="success">
                    <i className="fas fa-check-circle"></i>
                    {success}
                </div>
            )}

            <FileUpload onFileUpload={onFileUpload} disabled={loading} />
            
            <SampleProteins 
                samples={samples} 
                onSampleLoad={onSampleLoad} 
                disabled={loading}
            />
            
            {proteinData && (
                <>
                    <VisualizationControls
                        visualizationStyle={visualizationStyle}
                        setVisualizationStyle={setVisualizationStyle}
                        colorScheme={colorScheme}
                        setColorScheme={setColorScheme}
                    />
                    
                    <ProteinInfo proteinData={proteinData} />
                </>
            )}
        </div>
    );
}

// File Upload Component
function FileUpload({ onFileUpload, disabled }) {
    const fileInputRef = useRef();
    const [dragOver, setDragOver] = useState(false);

    const handleFileSelect = (file) => {
        if (file && (file.name.endsWith('.pdb') || file.name.endsWith('.ent'))) {
            onFileUpload(file);
        } else {
            alert('Please select a PDB file (.pdb or .ent)');
        }
    };

    const handleDrop = (e) => {
        e.preventDefault();
        setDragOver(false);
        
        const files = e.dataTransfer.files;
        if (files.length > 0) {
            handleFileSelect(files[0]);
        }
    };

    const handleDragOver = (e) => {
        e.preventDefault();
        setDragOver(true);
    };

    const handleDragLeave = (e) => {
        e.preventDefault();
        setDragOver(false);
    };

    const handleClick = () => {
        if (!disabled) {
            fileInputRef.current.click();
        }
    };

    const handleFileInputChange = (e) => {
        const file = e.target.files[0];
        if (file) {
            handleFileSelect(file);
        }
    };

    return (
        <div className="sidebar-section">
            <h3>
                <i className="fas fa-upload"></i>
                Upload PDB File
            </h3>
            
            <div
                className={`upload-area ${dragOver ? 'drag-over' : ''}`}
                onDrop={handleDrop}
                onDragOver={handleDragOver}
                onDragLeave={handleDragLeave}
                onClick={handleClick}
                style={{ opacity: disabled ? 0.6 : 1, cursor: disabled ? 'not-allowed' : 'pointer' }}
            >
                <i className="fas fa-cloud-upload-alt"></i>
                <p>Drag & drop your PDB file here</p>
                <p>or click to browse</p>
                <small>Supports .pdb and .ent files</small>
                
                <input
                    ref={fileInputRef}
                    type="file"
                    className="file-input"
                    accept=".pdb,.ent"
                    onChange={handleFileInputChange}
                    disabled={disabled}
                />
            </div>
        </div>
    );
}

// Sample Proteins Component
function SampleProteins({ samples, onSampleLoad, disabled }) {
    return (
        <div className="sidebar-section">
            <h3>
                <i className="fas fa-vial"></i>
                Sample Proteins
            </h3>
            
            <div className="sample-list">
                {samples.length === 0 ? (
                    <p style={{ color: '#718096', fontSize: '0.875rem', textAlign: 'center', padding: '1rem' }}>
                        No sample proteins available
                    </p>
                ) : (
                    samples.map((sample) => (
                        <div
                            key={sample.name}
                            className="sample-item"
                            onClick={() => !disabled && onSampleLoad(sample.name)}
                            style={{ 
                                opacity: disabled ? 0.6 : 1, 
                                cursor: disabled ? 'not-allowed' : 'pointer' 
                            }}
                        >
                            <h4>{sample.title}</h4>
                            <p>{sample.filename}</p>
                        </div>
                    ))
                )}
            </div>
        </div>
    );
}

// Visualization Controls Component
function VisualizationControls({ 
    visualizationStyle, 
    setVisualizationStyle,
    colorScheme,
    setColorScheme
}) {
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
        <div className="sidebar-section">
            <h3>
                <i className="fas fa-palette"></i>
                Visualization
            </h3>
            
            <div className="control-group">
                <label>Style</label>
                <div className="control-buttons">
                    {styles.map((style) => (
                        <button
                            key={style.value}
                            className={`btn btn-small ${visualizationStyle === style.value ? 'active' : 'btn-secondary'}`}
                            onClick={() => setVisualizationStyle(style.value)}
                        >
                            <i className={style.icon}></i>
                            {style.label}
                        </button>
                    ))}
                </div>
            </div>
            
            <div className="control-group">
                <label>Color Scheme</label>
                <select 
                    value={colorScheme} 
                    onChange={(e) => setColorScheme(e.target.value)}
                    style={{
                        width: '100%',
                        padding: '0.5rem',
                        border: '1px solid #e2e8f0',
                        borderRadius: '6px',
                        background: 'white'
                    }}
                >
                    {colors.map((color) => (
                        <option key={color.value} value={color.value}>
                            {color.label}
                        </option>
                    ))}
                </select>
            </div>
        </div>
    );
}

// Protein Info Component
function ProteinInfo({ proteinData }) {
    if (!proteinData) return null;

    const { header, stats } = proteinData;

    return (
        <div className="sidebar-section">
            <h3>
                <i className="fas fa-info-circle"></i>
                Protein Information
            </h3>
            
            <div className="protein-info">
                {header.title && (
                    <div style={{ marginBottom: '1rem' }}>
                        <h4 style={{ fontSize: '0.875rem', marginBottom: '0.5rem' }}>
                            {header.title}
                        </h4>
                    </div>
                )}
                
                <div className="info-grid">
                    <div className="info-item">
                        <span className="info-label">Atoms</span>
                        <span className="info-value">{stats.atom_count.toLocaleString()}</span>
                    </div>
                    
                    <div className="info-item">
                        <span className="info-label">Chains</span>
                        <span className="info-value">{stats.chain_count}</span>
                    </div>
                    
                    <div className="info-item">
                        <span className="info-label">Residues</span>
                        <span className="info-value">{stats.residue_count}</span>
                    </div>
                    
                    <div className="info-item">
                        <span className="info-label">Bonds</span>
                        <span className="info-value">{stats.bond_count.toLocaleString()}</span>
                    </div>
                    
                    {header.resolution && (
                        <div className="info-item">
                            <span className="info-label">Resolution</span>
                            <span className="info-value">{header.resolution} Å</span>
                        </div>
                    )}
                    
                    {header.experiment_type && (
                        <div className="info-item">
                            <span className="info-label">Method</span>
                            <span className="info-value">{header.experiment_type}</span>
                        </div>
                    )}
                    
                    {header.organism && (
                        <div className="info-item">
                            <span className="info-label">Organism</span>
                            <span className="info-value">{header.organism}</span>
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
}

// Protein Viewer Component
function ProteinViewer({ proteinData, visualizationStyle, colorScheme, loading }) {
    const viewerRef = useRef();
    const molViewerRef = useRef();
    const currentModelRef = useRef();

    // Initialize 3Dmol viewer
    useEffect(() => {
        if (viewerRef.current && !molViewerRef.current) {
            molViewerRef.current = $3Dmol.createViewer(viewerRef.current, {
                defaultcolors: $3Dmol.rasmolElementColors
            });
            
            // Set initial background
            molViewerRef.current.setBackgroundColor(0x1e3a8a);
        }
    }, []);

    // Load new protein data
    useEffect(() => {
        if (molViewerRef.current && proteinData) {
            loadProteinData();
        }
    }, [proteinData]);

    // Update visualization style when style or color changes
    useEffect(() => {
        if (molViewerRef.current && currentModelRef.current && proteinData) {
            applyVisualizationStyle(currentModelRef.current, visualizationStyle, colorScheme);
            molViewerRef.current.render();
        }
    }, [visualizationStyle, colorScheme]);

    const loadProteinData = () => {
        const viewer = molViewerRef.current;
        if (!viewer || !proteinData) return;

        // Clear existing models
        viewer.clear();
        currentModelRef.current = null;

        try {
            // Convert our protein data to PDB format for 3Dmol
            const pdbData = convertToPDBFormat(proteinData);
            
            // Add model to viewer
            const model = viewer.addModel(pdbData, 'pdb');
            currentModelRef.current = model;

            // Apply initial visualization style
            applyVisualizationStyle(model, visualizationStyle, colorScheme);

            // Center and zoom to fit
            viewer.zoomTo();
            viewer.render();
        } catch (error) {
            console.error('Error loading protein:', error);
        }
    };

    const convertToPDBFormat = (data) => {
        let pdbLines = [];
        
        // Add header information
        if (data.header.title) {
            pdbLines.push(`TITLE     ${data.header.title}`);
        }
        
        // Add atom records
        data.atoms.forEach(atom => {
            const line = [
                'ATOM  ',
                atom.serial.toString().padStart(5),
                '  ',
                atom.name.padEnd(4),
                atom.alt_loc || ' ',
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

    const applyVisualizationStyle = (model, style, color) => {
        const viewer = molViewerRef.current;
        if (!viewer) return;

        // Clear all existing styles on all models
        viewer.setStyle({}, {});

        // Resolve color scheme
        const scheme = (() => {
            switch (color) {
                case 'spectrum': return 'spectrum';
                case 'chain': return 'chain';
                case 'residue': return 'residue';
                case 'element': return 'default';
                default: return 'spectrum';
            }
        })();

        // Apply visualization style
        switch (style) {
            case 'cartoon':
                viewer.setStyle({}, {
                    cartoon: {
                        color: scheme === 'spectrum' ? 'spectrum' : scheme,
                        thickness: 0.4,
                        opacity: 0.9
                    }
                });
                break;
            case 'stick':
                viewer.setStyle({}, {
                    stick: {
                        colorscheme: scheme,
                        radius: 0.25
                    }
                });
                break;
            case 'sphere':
                viewer.setStyle({}, {
                    sphere: {
                        colorscheme: scheme,
                        scale: 0.8,
                        opacity: 0.9
                    }
                });
                break;
            case 'line':
                // CA trace to represent backbone cleanly
                viewer.setStyle({ atom: 'CA' }, {
                    line: {
                        colorscheme: scheme,
                        linewidth: 2
                    }
                });
                break;
            default:
                viewer.setStyle({}, {
                    cartoon: {
                        color: scheme === 'spectrum' ? 'spectrum' : scheme,
                        thickness: 0.4,
                        opacity: 0.9
                    }
                });
        }

        viewer.render();
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
            <div className="viewer-container">
                <div className="protein-viewer">
                    <div className="loading">
                        <div className="loading-spinner"></div>
                        <p>Loading protein structure...</p>
                    </div>
                </div>
            </div>
        );
    }

    return (
        <div className="viewer-container">
            <div 
                ref={viewerRef} 
                className="protein-viewer"
                style={{ 
                    width: '100%', 
                    height: '100%',
                    position: 'relative'
                }}
            >
                {!proteinData && (
                    <div className="loading">
                        <i className="fas fa-dna" style={{ fontSize: '3rem', marginBottom: '1rem' }}></i>
                        <p>Upload a PDB file or select a sample to begin</p>
                    </div>
                )}
            </div>
            
            {proteinData && (
                <div className="viewer-controls">
                    <button 
                        className="btn"
                        onClick={handleReset}
                        title="Reset View"
                    >
                        <i className="fas fa-home"></i>
                    </button>
                    
                    <button 
                        className="btn"
                        onClick={handleFullscreen}
                        title="Fullscreen"
                    >
                        <i className="fas fa-expand"></i>
                    </button>
                </div>
            )}
        </div>
    );
}

// Render the app
ReactDOM.render(<App />, document.getElementById('root'));
