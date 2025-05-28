# Repository for "A Competitive Disinhibitory Network for Robust Optic Flow Processing in Drosophila"

**Authors:**  
Mert Erginkaya $^{1,2}$, Tomás Cruz $^{1,3}$, Margarida Brotas $^{1,4}$, André Marques $^{1}$, Kathrin Steck $^{5}$, Aljoscha Nern $^{6}$, Filipa Torrão $^{1}$, Nélia Varela $^{1}$, Davi D. Bock $^{7}$, Michael Reiser $^{6}$, and M. Eugenia Chiappe $^{1,8,*}$

$^{1}$ Champalimaud Neuroscience Programme, Champalimaud Centre for the Unknown, Lisbon, Portugal  
$^{2}$ Neurobiology and Genetics, Theodor-Boveri-Institute, Biocenter, Julius-Maximilians-University of Würzburg, Würzburg, Germany  
$^{3}$ Friedrich Miescher Institute for Biomedical Research, and Biozentrum, Department of Cell Biology, University of Basel, Basel, Switzerland  
$^{4}$ CEDOC, iNOVA4Health, NOVA Medical School, Universidade Nova de Lisboa, Lisbon, Portugal  
$^{5}$ Center for Teacher Education and School Research, Universität Leipzig, Leipzig, Germany  
$^{6}$ Janelia Research Campus, Howard Hughes Medical Institute, Ashburn, VA  
$^{7}$ University of Vermont Larner College of Medicine, Burlington, VT  
$^{8}$ Lead contact

$^{*}$ Correspondence should be addressed to [eugenia.chiappe@neuro.fchampalimaud.org](mailto:eugenia.chiappe@neuro.fchampalimaud.org)

**Published in:**  
*Nature Neuroscience, 2025*\
**DOI:** [10.1038/s41593-025-01948-9](https://doi.org/10.1038/s41593-025-01948-9)

**Zenodo data repository:**\
[https://doi.org/10.5281/zenodo.14967806](https://doi.org/10.5281/zenodo.14967806)

---

## Overview
This repository contains the accompanying code and raw data associated with the publication “A Competitive Disinhibitory Network for Robust Optic Flow Processing in Drosophila.” 

---

## 📁 Repository Structure

```
| scripts/
| ├── matlab/                                  # Folder containing MATLAB scripts
| |    ├── PaperFiguresImaging.m                 # Main analysis and plotting script
| |    ├── PaperFiguresBehavior.m                # Behavior analysis and plotting script
| |    ├── Figure1b.m                            # Figure 1b plotting script
| |    └── helper-functions/                     # Helper functions for the scripts
| |
| └── r/                                       # Folder Containing .R scripts used for figures 
|      ├── 01_neuropil_analysis_outputs.Rmd      # R script for neuropil outputs analysis
|      ├── 02_neuropil_analysis_inputs.Rmd       # R script for neuropil inputs analysis 
|      ├── 03_proportion_diads_triads.Rmd        # R script for analysis of diads and triads
|      ├── 04_rootpoint_dist.Rmd                 # R script for rootpoint distance analysis
|      └── README.md                             # Acknowledgements and instructions for R scripts
|
| data/                                       # Folder containing data files located in Zenodo
| ├── Behavior/                                  # Behavior data for matlab scripts
| └── Preprocessed Data/                         # Preprocessed data for matlab scripts
|
| download_data_linux.sh                         # Script to download data on Linux/macOS
| download_data_windows.bat                      # Script to download data on Windows
| README.md
| License
```

---

## 🛠️ Instructions

### Step 1: Clone the Repository

You can clone the repository in two ways:

#### 1.1 Command Line

```bash
git clone https://github.com/your-username/your-publication-repo.git
cd your-publication-repo
```

#### 1.2 GitHub Desktop

1. Open GitHub Desktop.  
2. Click **File** → **Clone repository**.  
3. Select the repository from your list or paste the URL (above).  
4. Choose the local path where you want to clone it.

---

### Step 2: Download large Data

Large data files are hosted externally in Zenodo and can be downloaded using the provided download scripts for linux or windows based systems. Users can also download manually from the Zenodo repository and place the contents in the repository directory.

**Zenodo data repository:**\
[https://doi.org/10.5281/zenodo.14967806](https://doi.org/10.5281/zenodo.14967806) 

#### 2.1 Linux/macOS

1. Make sure you have `curl` or `wget` installed (for most Linux distributions, they are pre-installed).
2. Make sure you have `unzip` installed if you are downloading and extracting zipped folders.

   - **Debian/Ubuntu**:
     ```bash
     sudo apt-get update
     sudo apt-get install unzip -y
     ```
   - **macOS (Homebrew)**:
     ```bash
     brew install unzip
     ```

3. Navigate to the Repo directory and make the script executable and run it:
   ```bash
   chmod +x download_data_linux.sh
   ./download_data_linux.sh
   ```
   This will download and place the files into the appropriate folders (e.g. `data/Preprocessed Data/`).

#### 2.2 Windows

1. Double click the `download_data_windows.bat` press any key if prompted and wait for Download to finish

This will download and place the files into the appropriate folders (e.g. `data/Preprocessed Data/`).



---

### Step 3: Run the desired scripts

Once the data is downloaded and placed into the correct folders, you can open MATLAB or R scripts, and run the analysis. 

Note:
Be sure to add repository path to MATLAB path due to the use of helper functions.

---

## 🐝 License

This work is licensed under the [MIT License](License).

---

## 📧 Contact

For any questions or issues, please contact:  
M. Eugenia Chiappe: [eugenia.chiappe@neuro.fchampalimaud.org](mailto:eugenia.chiappe@neuro.fchampalimaud.org)

