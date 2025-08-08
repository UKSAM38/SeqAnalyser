import sys
import os
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                           QPushButton, QFileDialog, QTextEdit, QLabel,
                           QHBoxLayout, QListWidget, QProgressBar, QMessageBox, QScrollArea)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from sequence_utils import SequenceAnalyzer
from Bio import SeqIO

class AnalysisWorker(QThread):
    progress = pyqtSignal(int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)
    
    def __init__(self, analyzer, ab1_files, reference_seq, analyze_protein=False):
        super().__init__()
        self.analyzer = analyzer
        self.ab1_files = ab1_files
        self.reference_seq = reference_seq
        self.analyze_protein = analyze_protein
        
    def run(self):
        try:
            results = {'nucleotide': {}, 'protein': {}} if self.analyze_protein else {'nucleotide': {}}
            total_files = len(self.ab1_files)
            
            for i, ab1_file in enumerate(self.ab1_files):
                try:
                    # Read AB1 file
                    record = self.analyzer.read_ab1(ab1_file)
                    
                    # Trim sequence
                    trimmed_seq = self.analyzer.trim_sequence(record.seq)
                    
                    # Get filename without path and extension
                    filename = os.path.splitext(os.path.basename(ab1_file))[0]
                    
                    # Perform nucleotide alignment
                    nuc_alignment = self.analyzer.perform_alignment(trimmed_seq, self.reference_seq, filename, protein=False)
                    results['nucleotide'][filename] = nuc_alignment
                    
                    # Perform protein alignment if requested
                    if self.analyze_protein:
                        prot_alignment = self.analyzer.perform_alignment(trimmed_seq, self.reference_seq, filename, protein=True)
                        results['protein'][filename] = prot_alignment
                    
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
        self.setGeometry(100, 100, 1200, 800)
        
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
        
        # Create main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)
        main_layout.setSpacing(20)
        main_layout.setContentsMargins(20, 20, 20, 20)
        
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
        
        self.load_ab1_button = QPushButton("Load .ab1 Files")
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
        
        # Right panel (output)
        right_panel = QWidget()
        right_panel.setStyleSheet("""
            QWidget {
                background-color: #1E1B2E;
                border-radius: 8px;
            }
        """)
        right_layout = QVBoxLayout(right_panel)
        right_layout.setSpacing(10)
        right_layout.setContentsMargins(15, 15, 15, 15)
        
        results_title = QLabel("Analysis Results")
        results_title.setStyleSheet("font-size: 16px; color: #2D5BF5; margin-bottom: 10px;")
        right_layout.addWidget(results_title)
        
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
        right_layout.addWidget(scroll)
        
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
        right_layout.addWidget(self.progress_bar)
        
        # Add panels to main layout with proportions
        main_layout.addWidget(left_panel, 1)
        main_layout.addWidget(middle_panel, 2)
        main_layout.addWidget(right_panel, 3)
        
        # Initialize sequence analyzer and variables
        self.analyzer = SequenceAnalyzer()
        self.ab1_files = []
        self.reference_seq = None
        self.analysis_results = None
        
    def load_ab1_files(self):
        """Load AB1 files"""
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select AB1 Files",
            "",
            "AB1 Files (*.ab1)"
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
            "Sequence Files (*.txt *.fasta *.gb)"
        )
        
        if file:
            try:
                with open(file, 'r') as f:
                    self.reference_seq = f.read().strip()
                self.analyzer.set_reference(self.reference_seq)
                self.results_text.append("Reference sequence loaded")
                self.update_analyze_button()
            except Exception as e:
                self.results_text.append(f"Error loading reference: {str(e)}")
    
    def update_analyze_button(self):
        self.analyze_button.setEnabled(bool(self.ab1_files and self.reference_seq))
        
    def start_analysis(self, analyze_protein=False):
        if not hasattr(self, 'ab1_files') or not self.ab1_files:
            QMessageBox.warning(self, "Warning", "Please load .ab1 files first")
            return
            
        if not hasattr(self, 'reference_seq') or not self.reference_seq:
            QMessageBox.warning(self, "Warning", "Please load reference sequence first")
            return
            
        # Disable buttons during analysis
        self.load_ab1_button.setEnabled(False)
        self.load_ref_button.setEnabled(False)
        self.analyze_button.setEnabled(False)
        self.analyze_protein_button.setEnabled(False)
        self.export_pdf_button.setEnabled(False)
        
        # Clear previous results
        self.results_text.clear()
        
        # Create and start worker thread
        self.analyzer = SequenceAnalyzer()
        self.analyzer.set_reference(self.reference_seq)
        
        self.worker = AnalysisWorker(self.analyzer, self.ab1_files, self.reference_seq, analyze_protein)
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
        self.save_button.setEnabled(True)  # Enable save button
        
        # Store results for PDF export
        self.last_results = results
        
        # Display results
        for analysis_type, samples in results.items():
            self.results_text.append(f"\n=== {analysis_type.capitalize()} Analysis Results ===\n")
            for sample_name, result in samples.items():
                self.results_text.append(result)
                self.results_text.append("\n" + "="*50 + "\n")
        
        # Reset progress bar
        self.progress_bar.setVisible(False)
        
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

def main():
    app = QApplication(sys.argv)
    window = SeqAnalyseApp()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
