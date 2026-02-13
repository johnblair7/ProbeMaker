#!/usr/bin/env python3
"""
ProbeMaker Web Application
A web interface for generating complementary sequences to mRNA for given genes.
"""

from flask import Flask, render_template, request, send_file, jsonify
from flask_cors import CORS
import os
import tempfile
import time
from probe_maker import (
    ProbeMaker,
    GeneSequenceFetcher,
    apply_flex_handles,
    process_guide_rna_20mers,
    GUIDE_MODE_DEFAULT_LHS_25,
)

app = Flask(__name__)

# Enable CORS for all routes
CORS(app, origins=['*'])

# Global variables for session management
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'your-secret-key-here')
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size

@app.route('/')
def index():
    """Main page with gene input form."""
    return render_template('index.html')

@app.route('/generate', methods=['POST'])
def generate_probes():
    """Generate probe sequences: either from gene names (Probe mode) or from 20-base guide RNAs (Guide RNA mode)."""
    try:
        mode = request.form.get('mode', 'probe').strip().lower()
        flex_mode = request.form.get('flex_mode', 'v1').strip().lower()
        if flex_mode not in ('v1', 'v2'):
            flex_mode = 'v1'

        # --- Guide RNA mode: 20-base sequences only, no lookup ---
        if mode == 'guide_rna':
            guide_input = request.form.get('guide_sequences', '').strip()
            if not guide_input:
                return jsonify({'error': 'Please enter 20-base guide RNA sequences (one per line)'}), 400
            lines = [line.strip() for line in guide_input.split('\n') if line.strip()]
            if len(lines) > 500:
                return jsonify({'error': 'Maximum 500 guide sequences per request'}), 400
            try:
                results = process_guide_rna_20mers(lines, flex_mode=flex_mode, default_lhs_25=GUIDE_MODE_DEFAULT_LHS_25)
            except ValueError as e:
                return jsonify({'error': str(e)}), 400
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8') as temp_file:
                temp_file.write("LHS Probe (with handles)\tRHS Probe (with handles)\tGuide (20 bases)\n")
                temp_file.write("-" * 60 + "\t" + "-" * 60 + "\t" + "-" * 22 + "\n")
                for r in results:
                    temp_file.write(f"{r['lhs_with_handles']}\t{r['rhs_with_handles']}\t{r['guide_20']}\n")
                temp_file_path = temp_file.name
            return jsonify({
                'success': True,
                'message': f'Generated {len(results)} guide RNA probe rows (same LHS, unique RHS per guide)',
                'file_path': temp_file_path,
                'gene_count': len(results),
                'probe_count': len(results),
            })

        # --- Probe mode: gene names, lookup, probe design ---
        gene_input = request.form.get('gene_names', '').strip()
        num_pairs_raw = request.form.get('num_pairs', '').strip()
        try:
            num_pairs = int(num_pairs_raw) if num_pairs_raw else 3
        except ValueError:
            num_pairs = 3
        if num_pairs not in (2, 3):
            num_pairs = 3
        species = request.form.get('species', 'human').strip().lower()
        if species not in ('human', 'mouse'):
            species = 'human'
        run_blast = request.form.get('run_blast', '').strip().lower() in ('on', 'true', '1', 'yes')
        if not gene_input:
            return jsonify({'error': 'Please enter gene names'}), 400
        gene_names = [line.strip() for line in gene_input.split('\n') if line.strip()]
        if not gene_names:
            return jsonify({'error': 'No valid gene names found'}), 400
        if len(gene_names) > 100:
            return jsonify({'error': 'Maximum 100 genes allowed per request'}), 400
        probe_maker = ProbeMaker()
        sequence_fetcher = GeneSequenceFetcher()
        try:
            results = probe_maker.process_gene_names(gene_names, sequence_fetcher, num_pairs=num_pairs, species=species)
            if not results:
                return jsonify({'error': 'No valid probe pairs could be generated'}), 400
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8') as temp_file:
                temp_file.write("LHS Probe (25 bases)\tRHS Probe (25 bases)\tCombined Probe (50 bases)\tGene Name\n")
                temp_file.write("-" * 50 + "\t" + "-" * 50 + "\t" + "-" * 50 + "\t" + "-" * 20 + "\n")
                for result in results:
                    lhs_with_handles, rhs_with_handles = apply_flex_handles(
                        result['lhs_probe'], result['rhs_probe'], flex_mode
                    )
                    combined_probe = result['lhs_probe'] + result['rhs_probe']
                    gene_name = result.get('gene_name', 'Unknown')
                    temp_file.write(f"{lhs_with_handles}\t{rhs_with_handles}\t{combined_probe}\t{gene_name}\n")
                temp_file_path = temp_file.name
            if run_blast and results:
                try:
                    blast_report = probe_maker.generate_blast_report(results, species=species, timeout_seconds=30)
                    with open(temp_file_path, 'a', encoding='utf-8') as f:
                        f.write("\n\n")
                        f.write(blast_report)
                except Exception as blast_err:
                    with open(temp_file_path, 'a', encoding='utf-8') as f:
                        f.write(f"\n\nBLAST report skipped: {blast_err}\n")
            return jsonify({
                'success': True,
                'message': f'Successfully generated {len(results)} probe pairs for {len(gene_names)} genes',
                'file_path': temp_file_path,
                'gene_count': len(gene_names),
                'probe_count': len(results),
                'num_pairs_per_gene': num_pairs
            })
        finally:
            sequence_fetcher.close()
    except Exception as e:
        return jsonify({'error': f'Error generating probes: {str(e)}'}), 500

@app.route('/download/<path:filename>')
def download_file(filename):
    """Download the generated probe file."""
    try:
        # Security check - only allow .txt files from temp directory
        if not filename.endswith('.txt') or '..' in filename:
            return jsonify({'error': 'Invalid file'}), 400
        
        file_path = os.path.join(tempfile.gettempdir(), filename)
        
        if not os.path.exists(file_path):
            return jsonify({'error': 'File not found'}), 404
        
        # Send file and then delete it
        response = send_file(file_path, as_attachment=True, download_name='probe_sequences.txt')
        
        # Schedule file deletion after response is sent
        def cleanup():
            time.sleep(1)  # Wait for response to be sent
            try:
                os.unlink(file_path)
            except:
                pass
        
        # Start cleanup in background
        import threading
        cleanup_thread = threading.Thread(target=cleanup)
        cleanup_thread.daemon = True
        cleanup_thread.start()
        
        return response
        
    except Exception as e:
        return jsonify({'error': f'Error downloading file: {str(e)}'}), 500

@app.route('/health')
def health_check():
    """Health check endpoint for monitoring."""
    return jsonify({'status': 'healthy', 'service': 'ProbeMaker Web App'})

if __name__ == '__main__':
    # For development
    app.run(debug=True, host='0.0.0.0', port=8080)
else:
    # For production (e.g., gunicorn)
    pass
