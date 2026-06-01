
<h1>MYC Gene: Integrated Bioinformatics Analysis, Structural Validation, and Molecular Docking Pipeline</h1>

<div class="box">
<h2>Abstract</h2>
<p>The MYC proto-oncogene is a master regulator of cell growth, proliferation, apoptosis, and metabolism. Dysregulation of MYC is strongly associated with multiple human cancers, making it a critical target in oncology and drug discovery research. In this study, a complete in-silico bioinformatics pipeline was developed to analyze the MYC gene from sequence retrieval to molecular docking analysis. The workflow integrates nucleotide sequence extraction from NCBI, transcription into mRNA, translation into protein, sequence validation using BLAST, structural analysis using AlphaFold prediction, and molecular docking using GNINA. Additionally, ligand interaction analysis was performed using a PubChem-derived small molecule, and protein–ligand binding affinity was evaluated computationally. This project demonstrates how modern bioinformatics tools can be combined to simulate a complete drug discovery pipeline without wet-lab experiments.</p>
</div>

<div class="box" align="center">
<h2>Project WorkFlow</h2>

<div class="flowchart">

<div class="node">Retrieve MYC mRNA sequence from NCBI</div>
<div class="arrow">↓</div>

<div class="node">Extract coding DNA sequence (CDS)</div>
<div class="arrow">↓</div>

<div class="node">Transcribe DNA into mRNA</div>
<div class="arrow">↓</div>

<div class="node">Translate mRNA into protein sequence</div>
<div class="arrow">↓</div>

<div class="node">Validate sequence using BLAST</div>
<div class="arrow">↓</div>

<div class="node">Retrieve AlphaFold predicted structure</div>
<div class="arrow">↓</div>

<div class="node">Retrieve ligand from PubChem database</div>
<div class="arrow">↓</div>

<div class="node">Prepare ligand for docking</div>
<div class="arrow">↓</div>

<div class="node">Perform molecular docking using GNINA</div>
<div class="arrow">↓</div>

<div class="node">Analyze docking results</div>
<div class="arrow">↓</div>

<div class="node">Visualize protein–ligand interactions</div>

</div>
</div>

<div class="box">
<h2>Software and Tools Used</h2>

<table>
<tr><th>Tool</th><th>Description</th><th>Purpose</th></tr>
<tr><td>Biopython</td><td>Bioinformatics library</td><td>Sequence retrieval, transcription, translation, BLAST</td></tr>
<tr><td>NCBI Entrez</td><td>Database API</td><td>Fetch gene, mRNA, protein sequences</td></tr>
<tr><td>BLAST</td><td>Sequence alignment tool</td><td>Validate nucleotide/protein similarity</td></tr>
<tr><td>AlphaFold DB</td><td>AI structure prediction</td><td>Protein 3D structure prediction</td></tr>
<tr><td>Py3Dmol</td><td>Molecular visualization</td><td>3D structure visualization</td></tr>
<tr><td>NumPy</td><td>Numerical computing</td><td>Coordinate and docking center calculations</td></tr>
<tr><td>Pandas</td><td>Data analysis</td><td>Organizing BLAST and docking results</td></tr>
<tr><td>Matplotlib</td><td>Visualization</td><td>Docking score plots</td></tr>
<tr><td>PubChem</td><td>Chemical database</td><td>Ligand retrieval</td></tr>
<tr><td>Open Babel</td><td>Format conversion</td><td>Convert SDF → PDB → PDBQT</td></tr>
<tr><td>GNINA</td><td>AI docking engine</td><td>Protein–ligand docking</td></tr>
</table>

</div>

<div class="box">
<h2>Dataset Information</h2>

<table>
<tr><th>Component</th><th>Value</th></tr>
<tr><td>Gene</td><td>MYC (MYC proto-oncogene)</td></tr>
<tr><td>mRNA ID</td><td>NM_002467.6</td></tr>
<tr><td>Protein ID</td><td>NP_002458.2</td></tr>
<tr><td>AlphaFold Structure</td><td>AF-P01106-F1-model_v6</td></tr>
<tr><td>Ligand</td><td>PubChem Compound (CID used in code)</td></tr>
</table>

</div>

<div class="box">
<h2>Pipeline Description</h2>

<!-- 1 -->
<div class="box">
<h2>1. Environment Setup and Software Installation</h2>
<pre>
!pip install biopython requests pandas matplotlib py3Dmol -q
!apt-get update -y
!apt-get install openbabel -y -qq

!wget https://github.com/gnina/gnina/releases/download/v1.1/gnina -O gnina
!chmod +x gnina
!./gnina --help
</pre>
</div>

<P> This step installs and configures all required computational tools for the pipeline. Biopython, NumPy, Pandas, Matplotlib, and Py3Dmol are installed for bioinformatics analysis and visualization. Open Babel is installed for chemical format conversion. GNINA docking software is downloaded and made executable. Finally, the help command verifies successful installation and availability of docking parameters.</P>
<!-- 2 -->
<div class="box">
<h2>2. Importing Required Libraries</h2>
<pre>
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from IPython.display import display, HTML

import pandas as pd
import numpy as np
import requests
import py3Dmol
import matplotlib.pyplot as plt
import os
</pre>
</div>
<p>This section imports all required Python libraries for sequence analysis, structural biology, molecular docking, and visualization. Biopython modules handle sequence retrieval and structural parsing. NumPy and Pandas support numerical and tabular analysis. Py3Dmol enables 3D visualization of biomolecules, while Matplotlib is used for graphical representation of docking results.</p>
<!-- 3 -->
<div class="box">
<h2>3. Gene Definition and NCBI Sequence Retrieval</h2>
<pre>
Entrez.email = "subhanajamil3@gmail.com"

gene = "MYC"
ncbi_mrna_id = "NM_002467.6"
protein_ref_id = "NP_002458.2"

handle = Entrez.efetch(
    db="nucleotide",
    id="NM_002467.6",
    rettype="gb",
    retmode="text"
)

record = SeqIO.read(handle, "genbank")
handle.close()
</pre>
</div>
<p>The MYC gene is defined using its official NCBI accession numbers. The Entrez module is used to retrieve the full mRNA sequence in GenBank format. This format includes both nucleotide data and gene annotations. The sequence is parsed using SeqIO, and successful retrieval is confirmed through sequence metadata and length validation.</p>
<!-- 4 -->
<div class="box">
<h2>4. Coding Sequence (CDS) Extraction</h2>
<pre>
cds_seq = None

for feature in record.features:
    if feature.type == "CDS":
        cds_seq = feature.extract(record.seq)
        break

if cds_seq is None:
    raise ValueError("Coding region not found")
</pre>
</div>
<p>The GenBank file contains multiple annotated regions such as UTRs and CDS. This step isolates only the coding DNA sequence (CDS), which is responsible for protein synthesis. The CDS is extracted using feature-based parsing. Error handling ensures the pipeline stops if no coding region is found, maintaining biological validity.</p>
<!-- 5 -->
<div class="box">
<h2>5. mRNA Transcription</h2>
<pre>
mrna = cds_seq.transcribe()

print(mrna[:100])
</pre>
</div>
<p>The extracted CDS is transcribed into mRNA using Biopython’s transcription function. This simulates the biological process of gene expression where DNA is converted into RNA. During transcription, thymine bases are replaced by uracil. The resulting mRNA sequence maintains the same length as the CDS.</p>
<!-- 6 -->
<div class="box">
<h2>6. BLASTN Validation</h2>
<pre>
blast_result = NCBIWWW.qblast(
    "blastn",
    "nt",
    str(mrna),
    hitlist_size=5
)

with open("MYC_blast.xml", "w") as f:
    f.write(blast_result.read())
</pre>
</div>
<p>The transcribed mRNA sequence is validated using BLASTN against the NCBI nucleotide database. This step confirms biological authenticity by identifying homologous sequences. The top five hits are retrieved and stored in XML format for downstream analysis and comparison.</p>
<!-- 7 -->
<div class="box">
<h2>7. Protein Translation</h2>
<pre>
protein = mrna.translate(to_stop=True)

print(protein[:100])
</pre>
</div>
<p>The mRNA sequence is translated into a protein sequence using the standard genetic code. Each codon is converted into its corresponding amino acid. The process stops at the first stop codon, ensuring biologically relevant protein formation. The resulting MYC protein contains 454 amino acids.</p>
<!-- 8 -->
<div class="box">
<h2>8. Reference Protein Retrieval</h2>
<pre>
handle = Entrez.efetch(
    db="protein",
    id="NP_002458.2",
    rettype="fasta",
    retmode="text"
)

ref_protein = SeqIO.read(handle, "fasta").seq
handle.close()
</pre>
</div>
<p>The canonical MYC protein sequence is retrieved from the NCBI protein database. This reference sequence is used to validate the translated protein. Sequence comparison ensures correct reading frame and biological accuracy of translation.</p>
<!-- 9 -->
<div class="box">
<h2>9. AlphaFold Structure Retrieval</h2>
<pre>
url = "https://alphafold.ebi.ac.uk/files/AF-P01106-F1-model_v6.pdb"
r = requests.get(url)

with open("MYC_alphafold.pdb", "w") as f:
    f.write(r.text)
</pre>
</div>
<p>The 3D structure of MYC is downloaded from the AlphaFold Protein Structure Database. This structure is AI-predicted and represents a high-confidence model of the protein. The file is stored in PDB format for structural analysis and docking preparation.</p>
<!-- 10 -->
<div class="box">
<h2>10. AlphaFold Visualization</h2>
<pre>
view = py3Dmol.view(width=800, height=400)

with open("MYC_alphafold.pdb") as f:
    view.addModel(f.read(), "pdb")

view.setStyle({"cartoon": {"color": "spectrum"}})
view.zoomTo()
view.show()
</pre>
</div>
<p>The AlphaFold structure is visualized using Py3Dmol in a 3D interactive environment. The protein is displayed in cartoon representation to highlight secondary structural elements. A spectrum color scheme improves structural interpretation and clarity.</p>
 <img width="896" height="442" alt="image" src="https://github.com/user-attachments/assets/01112c09-513d-49eb-8541-8937a5246d3c" />
 <p>For 3D visualization Click this link</p> 
 https://colab.research.google.com/drive/1q-KWXehftshvs-PxKH6U2TFFzaTLMI5G#scrollTo=SgdUlcEgsunS&fullscreenOutput=true

<!-- 11 -->
<div class="box">
<h2>11. Ligand Retrieval (PubChem)</h2>
<pre>
cid = "9835049"
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"

r = requests.get(url)

if "V2000" in r.text or "M  END" in r.text:
    with open("MYC_ligand.sdf", "w") as f:
        f.write(r.text)
</pre>
</div>
<p>A small molecule ligand is retrieved from the PubChem database using its CID. The compound is downloaded in SDF format, which contains full structural information. File validation ensures correct chemical structure before proceeding to docking preparation.</p>
<!-- 12 -->
<div class="box">
<h2>12. Docking Center Calculation</h2>
<pre>
parser = PDBParser(QUIET=True)
structure = parser.get_structure("MYC", "MYC_alphafold.pdb")

coords = [atom.coord for atom in structure.get_atoms()]
coords = np.array(coords)

center = coords.mean(axis=0)
</pre>
</div>
<p>The docking center is calculated by extracting atomic coordinates from the AlphaFold structure. The mean coordinate of all atoms defines the geometric center of the protein. This point is used as the reference for defining the docking search box in GNINA</p>
<!-- 13 -->
<div class="box">
<h2>13. Molecular Docking using GNINA</h2>
<pre>
gnina_command = f"""./gnina \
-r MYC_alphafold.pdb \
-l MYC_ligand.pdbqt \
--center_x {center[0]:.3f} \
--center_y {center[1]:.3f} \
--center_z {center[2]:.3f} \
--size_x 40 \
--size_y 40 \
--size_z 40 \
--num_modes 10 \
-o MYC_docked.sdf
"""

!{gnina_command}
</pre>
</div>
<p>Molecular docking is performed using GNINA, a deep learning-enhanced docking engine. The AlphaFold structure is used as receptor and the prepared ligand as input. A cubic grid defines the search space, and 10 ligand poses are generated. Each pose is scored based on predicted binding affinity.</p>
<!-- 14 -->
<div class="box">
<h2>14. Docking Score Visualization</h2>
<pre>
import matplotlib.pyplot as plt

modes = list(range(1, 11))

scores = [
    -4.03, -3.82, -3.58, -3.56, -4.05,
    -3.58, -3.60, -4.16, -4.93, -3.61
]

plt.bar(modes, scores)
plt.gca().invert_yaxis()
plt.title("MYC GNINA Docking Scores (Affinity)")
plt.xlabel("Docking Mode")
plt.ylabel("Affinity (kcal/mol)")
plt.show()
</pre>
</div>
<p>Docking results are visualized using a bar plot to compare binding affinities across all poses. Lower (more negative) values indicate stronger binding. The plot reveals Pose 9 as the most stable binding conformation with the lowest energy score.</p>
<img width="637" height="514" alt="image" src="https://github.com/user-attachments/assets/15bbc192-9e13-43c9-b00f-3902507c5557" />

<!-- 15 -->
<div class="box">
<h2>15. Protein–Ligand Complex Visualization</h2>
<pre>
import py3Dmol

view = py3Dmol.view(width=800, height=400)

view.addModel(open("MYC_alphafold.pdb").read(), "pdb")
view.setStyle({"cartoon": {"color": "blue"}})

view.addModel(open("MYC_docked.sdf").read(), "sdf")
view.setStyle({"stick": {"colorscheme": "greenCarbon"}})

view.zoomTo()
view.show()
</pre>
</div>
<p>The final docked complex is visualized in 3D using Py3Dmol. The protein is displayed in cartoon form while the ligand is shown in stick representation. This visualization confirms ligand positioning within the predicted binding region of MYC, validating docking results.</p>
<img width="765" height="429" alt="image" src="https://github.com/user-attachments/assets/83b237f5-4c77-4630-b5ca-f1ecbbe0927a" />
<p>For 3D visualization Click this link</p>
https://colab.research.google.com/drive/1q-KWXehftshvs-PxKH6U2TFFzaTLMI5G#scrollTo=1UspbFL8y1aE&fullscreenOutput=true 
<div class="box">
<h2>Summary & Interpretation</h2>
<p>The pipeline demonstrates a complete in silico drug discovery workflow for the MYC protein. Although binding affinities indicate moderate interaction strength, the system provides a reproducible framework for future screening of stronger MYC inhibitors.</p>
</div>

<div class="box">
<h2>References</h2>
<ul>
<li>Biopython Documentation</li>
<li>NCBI</li>
<li>AlphaFold</li>
<li>PubChem</li>
<li>GNINA</li>
<li>Py3Dmol</li>
</ul>
</div>

</body>
</html>
