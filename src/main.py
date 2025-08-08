import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                           QPushButton, QFileDialog, QTextEdit, QLabel,
                           QHBoxLayout, QListWidget, QProgressBar, QMessageBox, QScrollArea,
                           QTabWidget, QComboBox, QSpinBox, QCheckBox, QGroupBox, QGridLayout,
                           QSplitter, QTreeWidget, QTreeWidgetItem, QTableWidget, QTableWidgetItem)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from sequence_utils import SequenceAnalyzer
from visualization import SequenceVisualizer
from phylogenetic import PhylogeneticAnalyzer
from primer_design import PrimerDesigner
from orf_analyzer import ORFAnalyzer
from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class AnalysisWorker(QThread):
    progress = pyqtSignal(int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    
    def __init__(self, analyzer, ab1_files, reference_seq, analyze_protein=False, advanced_options=None):
        super().__init__()
        self.analyzer = analyzer
        self.ab1_files = ab1_files
        self.reference_seq = reference_seq
        self.analyze_protein = analyze_protein
        self.advanced_options = advanced_options or {}
        
    def run(self):
        try:
            results = {'nucleotide': {}, 'protein': {}} if self.analyze_protein else {'nucleotide': {}}
            if self.advanced_options.get('include_statistics', False):
                results['statistics'] = {}
            if self.advanced_options.get('include_orf', False):
                results['orf_analysis'] = {}
            if self.advanced_options.get('include_restriction', False):
                results['restriction_sites'] = {}
                
            total_files = len(self.ab1_files)
            
            for i, ab1_file in enumerate(self.ab1_files):
                try:
                    # Read AB1 file
                    record = self.analyzer.read_multiple_formats(ab1_file)
                    
                    # Trim sequence
                    trimmed_seq = self.analyzer.trim_sequence(record.seq, record)
                    
                    # Get filename without path and extension
                    filename = os.path.splitext(os.path.basename(ab1_file))[0]
                    
                    # Perform nucleotide alignment
                    nuc_alignment = self.analyzer.perform_alignment(trimmed_seq, self.reference_seq, filename, protein=False)
                    results['nucleotide'][filename] = nuc_alignment
                    
                    # Perform protein alignment if requested
                    if self.analyze_protein:
                        prot_alignment = self.analyzer.perform_alignment(trimmed_seq, self.reference_seq, filename, protein=True)
                        results['protein'][filename] = prot_alignment
                        
                    # Additional analyses
                    if self.advanced_options.get('include_statistics', False):
                        stats = self.analyzer.calculate_sequence_statistics(trimmed_seq)
                        results['statistics'][filename] = stats
                        
                    if self.advanced_options.get('include_orf', False):
                        orf_analyzer = ORFAnalyzer()
                        summary, orfs = orf_analyzer.get_orf_summary(trimmed_seq)
                        results['orf_analysis'][filename] = {'summary': summary, 'orfs': orfs[:5]}  # Top 5 ORFs
                        
                    if self.advanced_options.get('include_restriction', False):
                        restriction_sites = self.analyzer.find_restriction_sites(trimmed_seq)
                        results['restriction_sites'][filename] = restriction_sites
                    
                    # Update progress
                    progress = int((i + 1) / total_files * 100)
                    self.progress.emit(progress)
                    
                except Exception as e:
                    self.error.emit(f"Error processing {ab1_file}: {str(e)}")
                    return
            
            self.finished.emit(results)
            
        except Exception as e:
            self.error.emit(str(e))

class SeqAnalyseApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SeqAnalyse - Bioinformatics Tool")
        self.setGeometry(100, 100, 1400, 900)
        
        # Define message box style for reuse
        self.message_box_style = """
            QMessageBox {
                background-color: #1E1B2E;
            }
            QMessageBox QLabel {
                color: #FFFFFF;
                font-size: 14px;
                padding: 10px;
            }
            QMessageBox QPushButton {
                background-color: #2D5BF5;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-size: 13px;
                min-width: 80px;
            }
            QMessageBox QPushButton:hover {
                background-color: #4171FF;
            }
        """
        
        self.setStyleSheet("""
            QMainWindow {
                background-color: #151421;
            }
            QWidget {
                color: #FFFFFF;
            }
            QPushButton {
                background-color: #2D5BF5;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                font-size: 14px;
                min-width: 120px;
                margin: 5px;
            }
            QPushButton:hover {
                background-color: #4171FF;
            }
            QPushButton:disabled {
                background-color: #2A2D3E;
                color: #6E7191;
            }
            QLabel {
                color: #FFFFFF;
                font-size: 14px;
                font-weight: bold;
            }
            QListWidget {
                background-color: #1E1B2E;
                border: 1px solid #2A2D3E;
                border-radius: 4px;
                padding: 5px;
                color: #FFFFFF;
            }
            QTabWidget::pane {
                border: 1px solid #2A2D3E;
                background-color: #1E1B2E;
            }
            QTabBar::tab {
                background-color: #2A2D3E;
                color: #FFFFFF;
                padding: 8px 16px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background-color: #2D5BF5;
            }
            QComboBox, QSpinBox {
                background-color: #2A2D3E;
                color: #FFFFFF;
                border: 1px solid #373B52;
                padding: 5px;
                border-radius: 4px;
            }
            QCheckBox {
                color: #FFFFFF;
            }
            QGroupBox {
                color: #FFFFFF;
                border: 1px solid #2A2D3E;
                border-radius: 4px;
                margin-top: 10px;
                padding-top: 10px;
            }
            QGroupBox::title {
                color: #2D5BF5;
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
            }
            QTreeWidget, QTableWidget {
                background-color: #1E1B2E;
                color: #FFFFFF;
                border: 1px solid #2A2D3E;
                border-radius: 4px;
            }
            QTreeWidget::item:selected, QTableWidget::item:selected {
                background-color: #2D5BF5;
            }
            QProgressBar {
                border: none;
                border-radius: 4px;
                background-color: #2A2D3E;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #2D5BF5;
                border-radius: 4px;
            }
            QScrollBar:vertical {
                border: none;
                background: #1E1B2E;
                width: 10px;
                border-radius: 5px;
            }
            QScrollBar::handle:vertical {
                background: #2A2D3E;
                border-radius: 5px;
            }
            QScrollBar::handle:vertical:hover {
                background: #373B52;
            }
            QTextEdit {
                background-color: #1E1B2E;
                color: #FFFFFF;
                border: none;
                border-radius: 4px;
                padding: 10px;
                selection-background-color: #2D5BF5;
            }
        """)
        
        self.setup_ui()
        self.setup_analyzers()
        
    def setup_analyzers(self):
        """Initialize all analyzers"""
        self.analyzer = SequenceAnalyzer()
        self.visualizer = SequenceVisualizer()
        self.phylo_analyzer = PhylogeneticAnalyzer()
        self.primer_designer = PrimerDesigner()
        self.orf_analyzer = ORFAnalyzer()
        self.ab1_files = []
        self.reference_seq = None
        self.analysis_results = None
        
    def setup_ui(self):
        """Setup the user interface"""
        # Create main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        
        # Create main splitter
        main_splitter = QSplitter(Qt.Orientation.Horizontal)
        main_layout = QHBoxLayout(main_widget)
        main_layout.addWidget(main_splitter)
        
        # Left panel (buttons)
        left_panel = QWidget()
        left_panel.setFixedWidth(200)
        left_panel.setStyleSheet("""
            QWidget {
                background-color: #1E1B2E;
                border-radius: 8px;
                padding: 10px;
            }
        """)
        left_layout = QVBoxLayout(left_panel)
        left_layout.setSpacing(10)
        left_layout.setContentsMargins(15, 15, 15, 15)
        
        # Add title to left panel
        title_label = QLabel("Controls")
        title_label.setStyleSheet("""
            font-size: 16px;
            color: #2D5BF5;
            padding: 5px;
            margin-bottom: 10px;
        """)
        left_layout.addWidget(title_label)
        
        # Create buttons with consistent styling
        button_style = """
            QPushButton {
                background-color: #2D5BF5;
                color: white;
                border: none;
                padding: 10px;
                border-radius: 4px;
                font-size: 13px;
                min-height: 35px;
            }
            QPushButton:hover {
                background-color: #4171FF;
            }
            QPushButton:disabled {
                background-color: #2A2D3E;
                color: #6E7191;
            }
        """
        
        self.load_ab1_button = QPushButton("Load Sequence Files")
        self.load_ab1_button.setStyleSheet(button_style)
        self.load_ab1_button.clicked.connect(self.load_ab1_files)
        left_layout.addWidget(self.load_ab1_button)
        
        self.load_ref_button = QPushButton("Load Reference")
        self.load_ref_button.setStyleSheet(button_style)
        self.load_ref_button.clicked.connect(self.load_reference)
        left_layout.addWidget(self.load_ref_button)
        
        # Add spacing between load and analyze buttons
        left_layout.addSpacing(10)
        
        self.analyze_button = QPushButton("Analyze Nucleotide")
        self.analyze_button.setStyleSheet(button_style)
        self.analyze_button.clicked.connect(lambda: self.start_analysis(False))
        left_layout.addWidget(self.analyze_button)
        
        self.analyze_protein_button = QPushButton("Analyze Protein")
        self.analyze_protein_button.setStyleSheet(button_style)
        self.analyze_protein_button.clicked.connect(lambda: self.start_analysis(True))
        left_layout.addWidget(self.analyze_protein_button)
        
        # Advanced analysis buttons
        self.phylo_button = QPushButton("Phylogenetic Tree")
        self.phylo_button.setStyleSheet(button_style)
        self.phylo_button.clicked.connect(self.create_phylogenetic_tree)
        self.phylo_button.setEnabled(False)
        left_layout.addWidget(self.phylo_button)
        
        self.primer_button = QPushButton("Design Primers")
        self.primer_button.setStyleSheet(button_style)
        self.primer_button.clicked.connect(self.design_primers)
        self.primer_button.setEnabled(False)
        left_layout.addWidget(self.primer_button)
        
        # Add spacing between analyze and export buttons
        left_layout.addSpacing(10)
        
        self.export_pdf_button = QPushButton("Export to PDF")
        self.export_pdf_button.setStyleSheet(button_style)
        self.export_pdf_button.clicked.connect(self.export_to_pdf)
        self.export_pdf_button.setEnabled(False)
        left_layout.addWidget(self.export_pdf_button)
        
        self.save_button = QPushButton("Save to Word")
        self.save_button.setStyleSheet(button_style)
        self.save_button.clicked.connect(self.save_results)
        self.save_button.setEnabled(False)
        left_layout.addWidget(self.save_button)
        
        self.export_html_button = QPushButton("Export to HTML")
        self.export_html_button.setStyleSheet(button_style)
        self.export_html_button.clicked.connect(self.export_to_html)
        self.export_html_button.setEnabled(False)
        left_layout.addWidget(self.export_html_button)
        
        left_layout.addStretch()
        
        # Middle panel (file list)
        middle_panel = QWidget()
        middle_panel.setMinimumWidth(300)
        middle_panel.setStyleSheet("""
            QWidget {
                background-color: #1E1B2E;
                border-radius: 8px;
            }
        """)
        middle_layout = QVBoxLayout(middle_panel)
        middle_layout.setSpacing(10)
        middle_layout.setContentsMargins(15, 15, 15, 15)
        
        files_title = QLabel("Loaded Files")
        files_title.setStyleSheet("""
            font-size: 16px;
            color: #2D5BF5;
            padding: 5px;
            margin-bottom: 10px;
        """)
        middle_layout.addWidget(files_title)
        
        self.file_list = QListWidget()
        self.file_list.setStyleSheet("""
            QListWidget {
                background-color: #1E1B2E;
                border: 1px solid #2A2D3E;
                border-radius: 4px;
                padding: 5px;
                font-size: 13px;
                color: #FFFFFF;
            }
            QListWidget::item {
                padding: 8px;
                border-bottom: 1px solid #2A2D3E;
            }
            QListWidget::item:selected {
                background-color: #2D5BF5;
                color: #FFFFFF;
            }
            QListWidget::item:hover {
                background-color: #2A2D3E;
            }
        """)
        middle_layout.addWidget(self.file_list)
        
        # Advanced options
        options_group = QGroupBox("Analysis Options")
        options_layout = QVBoxLayout(options_group)
        
        self.include_stats_cb = QCheckBox("Include sequence statistics")
        self.include_orf_cb = QCheckBox("Include ORF analysis")
        self.include_restriction_cb = QCheckBox("Include restriction sites")
        self.quality_threshold_spin = QSpinBox()
        self.quality_threshold_spin.setRange(10, 40)
        self.quality_threshold_spin.setValue(20)
        
        options_layout.addWidget(self.include_stats_cb)
        options_layout.addWidget(self.include_orf_cb)
        options_layout.addWidget(self.include_restriction_cb)
        options_layout.addWidget(QLabel("Quality threshold:"))
        options_layout.addWidget(self.quality_threshold_spin)
        
        middle_layout.addWidget(options_group)
        
        # Right panel (output)
        right_panel = QTabWidget()
        right_panel.setStyleSheet("""
            QTabWidget {
                background-color: #1E1B2E;
                border-radius: 8px;
            }
        """)
        
        # Results tab
        results_tab = QWidget()
        results_layout = QVBoxLayout(results_tab)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setStyleSheet("""
            QScrollArea {
                border: none;
            }
            QScrollBar:vertical {
                border: none;
                background: #1E1B2E;
                width: 10px;
                border-radius: 5px;
            }
            QScrollBar::handle:vertical {
                background: #2A2D3E;
                border-radius: 5px;
            }
            QScrollBar::handle:vertical:hover {
                background: #373B52;
            }
        """)
        results_layout.addWidget(scroll)
        
        results_widget = QWidget()
        scroll.setWidget(results_widget)
        self.results_layout = QVBoxLayout(results_widget)
        
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        font = QFont("Courier New", 10)
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setFixedPitch(True)
        self.results_text.setFont(font)
        self.results_text.setStyleSheet("""
            QTextEdit {
                background-color: #1E1B2E;
                color: #FFFFFF;
                border: none;
                border-radius: 4px;
                padding: 10px;
                selection-background-color: #2D5BF5;
            }
        """)
        self.results_layout.addWidget(self.results_text)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: none;
                border-radius: 4px;
                background-color: #2A2D3E;
                height: 8px;
                text-align: center;
            }
            QProgressBar::chunk {
                background-color: #2D5BF5;
                border-radius: 4px;
            }
        """)
        results_layout.addWidget(self.progress_bar)
        
        # Visualization tab
        viz_tab = QWidget()
        viz_layout = QVBoxLayout(viz_tab)
        viz_layout.addWidget(self.visualizer.get_canvas())
        
        # Phylogenetic tab
        phylo_tab = QWidget()
        phylo_layout = QVBoxLayout(phylo_tab)
        phylo_layout.addWidget(self.phylo_analyzer.get_canvas())
        
        # Statistics tab
        stats_tab = QWidget()
        stats_layout = QVBoxLayout(stats_tab)
        self.stats_tree = QTreeWidget()
        self.stats_tree.setHeaderLabels(["Property", "Value"])
        stats_layout.addWidget(self.stats_tree)
        
        # Add tabs
        right_panel.addTab(results_tab, "Results")
        right_panel.addTab(viz_tab, "Visualization")
        right_panel.addTab(phylo_tab, "Phylogenetic")
        right_panel.addTab(stats_tab, "Statistics")
        
        # Add panels to main layout with proportions
        main_splitter.addWidget(left_panel)
        main_splitter.addWidget(middle_panel)
        main_splitter.addWidget(right_panel)
        main_splitter.setSizes([200, 300, 900])
        
    def load_ab1_files(self):
        """Load sequence files"""
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select Sequence Files",
            "",
            "Sequence Files (*.ab1 *.fasta *.fa *.fastq *.fq *.gb *.gbk *.embl *.txt);;AB1 Files (*.ab1);;FASTA Files (*.fasta *.fa);;FASTQ Files (*.fastq *.fq);;GenBank Files (*.gb *.gbk);;EMBL Files (*.embl);;Text Files (*.txt)"
        )
        
        if files:
            self.ab1_files = files
            self.file_list.clear()
            for file in files:
                self.file_list.addItem(os.path.basename(file))
            self.update_analyze_button()
            
    def load_reference(self):
        """Load reference sequence file"""
        file, _ = QFileDialog.getOpenFileName(
            self,
            "Select Reference File",
            "",
            "Sequence Files (*.txt *.fasta *.fa *.gb *.gbk *.embl);;Text Files (*.txt);;FASTA Files (*.fasta *.fa);;GenBank Files (*.gb *.gbk);;EMBL Files (*.embl)"
        )
        
        if file:
            try:
                # Try to read as sequence file first
                try:
                    record = self.analyzer.read_multiple_formats(file)
                    self.reference_seq = str(record.seq)
                except:
                    # Fallback to plain text
                    with open(file, 'r') as f:
                        self.reference_seq = f.read().strip()
                        
                self.analyzer.set_reference(self.reference_seq)
                self.results_text.append("Reference sequence loaded")
                self.update_analyze_button()
            except Exception as e:
                self.results_text.append(f"Error loading reference: {str(e)}")
    
    def update_analyze_button(self):
        self.analyze_button.setEnabled(bool(self.ab1_files and self.reference_seq))
        self.analyze_protein_button.setEnabled(bool(self.ab1_files and self.reference_seq))
        
    def start_analysis(self, analyze_protein=False):
        if not self.ab1_files:
            QMessageBox.warning(self, "Warning", "Please load sequence files first")
            return
            
        if not self.reference_seq:
            QMessageBox.warning(self, "Warning", "Please load reference sequence first")
            return
            
        # Disable buttons during analysis
        self.load_ab1_button.setEnabled(False)
        self.load_ref_button.setEnabled(False)
        self.analyze_button.setEnabled(False)
        self.analyze_protein_button.setEnabled(False)
        self.export_pdf_button.setEnabled(False)
        self.phylo_button.setEnabled(False)
        self.primer_button.setEnabled(False)
        
        # Clear previous results
        self.results_text.clear()
        
        # Get advanced options
        advanced_options = {
            'include_statistics': self.include_stats_cb.isChecked(),
            'include_orf': self.include_orf_cb.isChecked(),
            'include_restriction': self.include_restriction_cb.isChecked(),
            'quality_threshold': self.quality_threshold_spin.value()
        }
        
        # Create and start worker thread
        self.analyzer = SequenceAnalyzer()
        self.analyzer.quality_threshold = advanced_options['quality_threshold']
        self.analyzer.set_reference(self.reference_seq)
        
        self.worker = AnalysisWorker(self.analyzer, self.ab1_files, self.reference_seq, analyze_protein, advanced_options)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.analysis_complete)
        self.worker.error.connect(self.analysis_error)
        self.worker.start()
        
    def update_progress(self, value):
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(value)
        
    def analysis_complete(self, results):
        # Enable buttons
        self.load_ab1_button.setEnabled(True)
        self.load_ref_button.setEnabled(True)
        self.analyze_button.setEnabled(True)
        self.analyze_protein_button.setEnabled(True)
        self.export_pdf_button.setEnabled(True)
        self.export_html_button.setEnabled(True)
        self.save_button.setEnabled(True)  # Enable save button
        self.phylo_button.setEnabled(True)
        self.primer_button.setEnabled(True)
        
        # Store results for PDF export
        self.last_results = results
        
        # Display results
        for analysis_type, samples in results.items():
            if analysis_type in ['nucleotide', 'protein']:
                self.results_text.append(f"\n=== {analysis_type.capitalize()} Analysis Results ===\n")
                for sample_name, result in samples.items():
                    self.results_text.append(result)
                    self.results_text.append("\n" + "="*50 + "\n")
                    
        # Update statistics tree
        self.update_statistics_tree(results)
        
        # Update visualizations
        if self.ab1_files:
            try:
                record = self.analyzer.read_multiple_formats(self.ab1_files[0])
                self.visualizer.plot_quality_scores(record)
            except:
                pass
        
        # Reset progress bar
        self.progress_bar.setVisible(False)
        
    def update_statistics_tree(self, results):
        """Update the statistics tree widget"""
        self.stats_tree.clear()
        
        if 'statistics' in results:
            for sample_name, stats in results['statistics'].items():
                sample_item = QTreeWidgetItem(self.stats_tree, [sample_name, ""])
                
                # Basic stats
                basic_item = QTreeWidgetItem(sample_item, ["Basic Statistics", ""])
                QTreeWidgetItem(basic_item, ["Length", str(stats['length'])])
                QTreeWidgetItem(basic_item, ["GC Content", f"{stats['gc_content']:.2f}%"])
                QTreeWidgetItem(basic_item, ["Molecular Weight", f"{stats['molecular_weight']:.2f} Da"])
                
                # Base composition
                comp_item = QTreeWidgetItem(sample_item, ["Base Composition", ""])
                for base, count in stats['base_counts'].items():
                    percentage = stats['base_percentages'][base]
                    QTreeWidgetItem(comp_item, [base, f"{count} ({percentage:.1f}%)"])
                    
                sample_item.setExpanded(True)
                basic_item.setExpanded(True)
                comp_item.setExpanded(True)
        
    def analysis_error(self, error_msg):
        msg = QMessageBox()
        msg.setStyleSheet(self.message_box_style)
        msg.setIcon(QMessageBox.Icon.Critical)
        msg.setText(f"Analysis failed: {error_msg}")
        msg.setWindowTitle("Error")
        msg.exec()
        self.analyze_button.setEnabled(True)
        self.analyze_protein_button.setEnabled(True)
        self.progress_bar.setVisible(False)
        
    def create_phylogenetic_tree(self):
        """Create phylogenetic tree from analyzed sequences"""
        if not hasattr(self, 'last_results') or 'nucleotide' not in self.last_results:
            QMessageBox.warning(self, "Warning", "No nucleotide analysis results available")
            return
            
        try:
            # Extract sequences from results
            sequences = {}
            for sample_name in self.last_results['nucleotide'].keys():
                # Get the trimmed sequence for each sample
                record = self.analyzer.read_multiple_formats(
                    next(f for f in self.ab1_files if os.path.splitext(os.path.basename(f))[0] == sample_name)
                )
                trimmed_seq = self.analyzer.trim_sequence(record.seq, record)
                sequences[sample_name] = trimmed_seq
                
            # Add reference sequence
            sequences['Reference'] = self.reference_seq
            
            # Create phylogenetic tree
            tree, distance_matrix = self.phylo_analyzer.analyze_sequences(sequences)
            
            QMessageBox.information(self, "Success", "Phylogenetic tree created successfully!")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create phylogenetic tree: {str(e)}")
            
    def design_primers(self):
        """Design PCR primers"""
        if not self.reference_seq:
            QMessageBox.warning(self, "Warning", "Please load reference sequence first")
            return
            
        try:
            primers = self.primer_designer.design_primers(self.reference_seq)
            primer_pairs = self.primer_designer.find_primer_pairs(primers)
            
            # Display primer results
            primer_text = "=== Primer Design Results ===\n\n"
            
            if primer_pairs:
                primer_text += "Top Primer Pairs:\n\n"
                for i, pair in enumerate(primer_pairs, 1):
                    primer_text += f"Pair {i}:\n"
                    primer_text += f"  Forward: {pair['forward']['sequence']}\n"
                    primer_text += f"    Tm: {pair['forward']['tm']:.1f}°C, GC: {pair['forward']['gc_content']:.1f}%\n"
                    primer_text += f"  Reverse: {pair['reverse']['sequence']}\n"
                    primer_text += f"    Tm: {pair['reverse']['tm']:.1f}°C, GC: {pair['reverse']['gc_content']:.1f}%\n"
                    primer_text += f"  Product size: {pair['product_size']} bp\n"
                    primer_text += f"  Quality score: {pair['quality_score']:.3f}\n\n"
            else:
                primer_text += "No suitable primer pairs found.\n\n"
                
            if primers:
                primer_text += "Individual Primers:\n\n"
                for i, primer in enumerate(primers[:10], 1):
                    primer_text += f"{primer['type'].capitalize()} {i}: {primer['sequence']}\n"
                    primer_text += f"  Tm: {primer['tm']:.1f}°C, GC: {primer['gc_content']:.1f}%, Length: {primer['length']} bp\n\n"
                    
            self.results_text.append(primer_text)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Primer design failed: {str(e)}")
        
    def save_results(self):
        """Save results to Word document"""
        if not hasattr(self, 'last_results'):
            msg = QMessageBox()
            msg.setStyleSheet(self.message_box_style)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText("No analysis results to save")
            msg.setWindowTitle("Warning")
            msg.exec()
            return
            
        try:
            file_name, _ = QFileDialog.getSaveFileName(
                self,
                "Save Results",
                "",
                "Word Files (*.docx)"
            )
            if file_name:
                if not file_name.endswith('.docx'):
                    file_name += '.docx'
                self.analyzer.save_to_word(self.last_results, file_name)
                self.results_text.append(f"\nResults saved to {file_name}")
        except Exception as e:
            msg = QMessageBox()
            msg.setStyleSheet(self.message_box_style)
            msg.setIcon(QMessageBox.Icon.Critical)
            msg.setText(f"Error saving to Word: {str(e)}")
            msg.setWindowTitle("Error")
            msg.exec()
            
    def export_to_pdf(self):
        if not hasattr(self, 'last_results'):
            msg = QMessageBox()
            msg.setStyleSheet(self.message_box_style)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText("No analysis results to export!")
            msg.setWindowTitle("Export Error")
            msg.exec()
            return
            
        # Get save file name
        file_name, _ = QFileDialog.getSaveFileName(self, "Save PDF", "", "PDF Files (*.pdf)")
        if not file_name:
            return
            
        if not file_name.endswith('.pdf'):
            file_name += '.pdf'
            
        try:
            self.analyzer.export_to_pdf(self.last_results, file_name)
            msg = QMessageBox()
            msg.setStyleSheet(self.message_box_style)
            msg.setIcon(QMessageBox.Icon.Information)
            msg.setText(f"Results exported to {file_name}")
            msg.setWindowTitle("Success")
            msg.exec()
        except Exception as e:
            msg = QMessageBox()
            msg.setStyleSheet(self.message_box_style)
            msg.setIcon(QMessageBox.Icon.Critical)
            msg.setText(f"Failed to export PDF: {str(e)}")
            msg.setWindowTitle("Error")
            msg.exec()
            
    def export_to_html(self):
        """Export results to HTML"""
        if not hasattr(self, 'last_results'):
            QMessageBox.warning(self, "Warning", "No analysis results to export!")
            return
            
        file_name, _ = QFileDialog.getSaveFileName(self, "Save HTML", "", "HTML Files (*.html)")
        if not file_name:
            return
            
        if not file_name.endswith('.html'):
            file_name += '.html'
            
        try:
            self.analyzer.export_advanced_results(self.last_results, file_name, 'html')
            QMessageBox.information(self, "Success", f"Results exported to {file_name}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export HTML: {str(e)}")

def main():
    app = QApplication(sys.argv)
    window = SeqAnalyseApp()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
