from flask import (Blueprint, redirect, render_template, request)
from flaskr.auth import login_required

import os

import rosalind


bp = Blueprint('blog', __name__)

bp.config = {}
bp.config["FILE_UPLOADS"] = "C:/Users/enesd/bioapp/Uploads"
bp.config["FILE_DOWNLOADS"] = "C:/Users/enesd/bioapp/Uploads"
bp.config["ALLOWED_FILE_EXTENSIONS"] = ["FASTA", "TXT"]



def allowed_file(filename):

    if not "." in filename:
        return False

    ext = filename.rsplit(".", 1)[1]

    if ext.upper() in bp.config["ALLOWED_FILE_EXTENSIONS"]:
        return True
    else:
        return False

@bp.route('/', methods = ['GET', 'POST'])
@login_required
def upload():
    fas={}
    func= ['DNAcount', 'DNAtoRNA', 'reverseComplement', 'GCContent', 'countMutations', 'RNAtoProtein', 'DNAmotif', 'consensus']
    if request.method == 'POST':
        
        fasta = request.form['fasta']
        
        if fasta  is not None:
            
            with open(bp.config['FILE_DOWNLOADS']+"/fasta.fasta", "w") as fasta_file:
                fasta_file.write(fasta)
                
            fas = rosalind.readFasta("C:/Users/enesd/Desktop/İMÜ/Kodlama/bioapp/Uploads/fasta.fasta")
            
                
        if request.files:
            
            files = request.files.getlist("file")
            
            for file in files:
                print(file.filename)
                if file.filename == "":
                    print("Dosya ismi yok.")
                    #return redirect(request.url)
                
                if allowed_file(file.filename):
                    
                    file.save(os.path.join(bp.config['FILE_UPLOADS'], 'fasta.fasta'))
                    
                    print("Dosya yüklemesi başarılı.")
                    fas = rosalind.readFasta("C:/Users/enesd/Desktop/İMÜ/Kodlama/bioapp/Uploads/fasta.fasta")
                    
                    #return redirect(request.url)
                
                else:
                    print("Bu dosya türü uygun değil.")
                    #return redirect(request.url)
        return redirect("/apply")

    
    return render_template("blog/upload.html", variable=fas, functions=func)

@bp.route('/apply', methods = ['GET', 'POST'])
@login_required
def apply():
    func= ['DNAcount', 'DNAtoRNA', 'reverseComplement', 'GCContent', 'countMutations', 'RNAtoProtein', 'DNAmotif', 'consensus']
    fas = rosalind.readFasta("C:/Users/enesd/Desktop/İMÜ/Kodlama/bioapp/Uploads/fasta.fasta")
    length = len(fas)
    result=0
    liste=[]
    for i in range(0, length):
        liste.append(f'seq{i}')
        
    if request.method == 'POST':
        if request.form.get('function') == 'DNAcount':
            result = rosalind.DNAcount(fas[request.form.getlist('seq')[0]])
        if request.form.get('function') == 'DNAtoRNA':
            result = rosalind.DNAtoRNA(fas[request.form.getlist('seq')[0]])
        if request.form.get('function') == 'reverseComplement':
            result = rosalind.reverseComplement(fas[request.form.getlist('seq')[0]])
        if request.form.get('function') == 'GCContent':
            result = rosalind.GCContent(fas[request.form.getlist('seq')[0]])
        if request.form.get('function') == 'countMutations':
            result = rosalind.countMutations(fas[request.form.getlist('seq')[0]], fas[request.form.getlist('seq')[1]])
        if request.form.get('function') == 'RNAtoProtein':
            result = rosalind.RNAtoProtein(fas[request.form.getlist('seq')[0]])
        if request.form.get('function') == 'DNAmotif':
            result = rosalind.DNAmotif(fas[request.form.getlist('seq')[0]], fas[request.form.getlist('seq')[1]])
        if request.form.get('function') == 'consensus':
            result = rosalind.consensus([fas[i] for i in request.form.getlist('seq')])
                        
    
    return render_template("blog/apply.html", variable=fas, functions=func, sequences=liste, sonuc=result)
  