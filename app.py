#!/usr/bin/env python3
"""
ProbeMaker Web Application
A web interface for generating complementary sequences to mRNA for given genes.
"""

from flask import Flask, render_template, request, send_file, jsonify
import os
import tempfile
import time
from probe_maker import ProbeMaker, GeneSequenceFetcher

app = Flask(__name__)

# Global variables for session management
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'your-secret-key-here')
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size

@app.route('/')
def index():
    """Main page with gene input form."""
    return render_template('index.html')

@app.route('/generate', methods=['POST'])
def generate_probes():
    """Generate probe sequences for the submitted gene names."""
    try:
        # Get gene names from form
        gene_input = request.form.get('gene_names', '').strip()
        
        if not gene_input:
            return jsonify({'error': 'Please enter gene names'}), 400
        
        # Parse gene names (one per line)
        gene_names = [line.strip() for line in gene_input.split('\n') if line.strip()]
        
        if not gene_names:
            return jsonify({'error': 'No valid gene names found'}), 400
        
        if len(gene_names) > 100:  # Limit to prevent abuse
            return jsonify({'error': 'Maximum 100 genes allowed per request'}), 400
        
        # Initialize ProbeMaker and GeneSequenceFetcher
        probe_maker = ProbeMaker()
        sequence_fetcher = GeneSequenceFetcher()
        
        try:
            # Process gene names and generate probe pairs
            results = probe_maker.process_gene_names(gene_names, sequence_fetcher)
            
            if not results:
                return jsonify({'error': 'No valid probe pairs could be generated'}), 400
            
            # Create temporary file for download
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8') as temp_file:
                # Write header
                temp_file.write("LHS Probe (25 bases)\tRHS Probe (25 bases)\tGene Name\n")
                temp_file.write("-" * 50 + "\t" + "-" * 50 + "\t" + "-" * 20 + "\n")
                
                # Write probe data
                for result in results:
                    lhs_probe = result['lhs_probe']
                    rhs_probe = result['rhs_probe']
                    gene_name = result.get('gene_name', 'Unknown')
                    temp_file.write(f"{lhs_probe}\t{rhs_probe}\t{gene_name}\n")
                
                temp_file_path = temp_file.name
            
            # Return success response with file info
            return jsonify({
                'success': True,
                'message': f'Successfully generated {len(results)} probe pairs for {len(gene_names)} genes',
                'file_path': temp_file_path,
                'gene_count': len(gene_names),
                'probe_count': len(results)
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
