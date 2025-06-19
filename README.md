# Project_Bioinformatics

Repository created as part of a Bioinformatics master's degree project "Development of a User-Friendly Graphical Platform for Molecular Characterization of Macroscopic Parasites in Fish".

Created by [João Faria](https://github.com/JohnnyFarians24) (PG55700).

---

# The PARAFISH Platform

PARAFISH is an interactive web-based platform built with Dash for the streamlined molecular analysis of sequencing data. It provides a step-by-step graphical workflow, guiding the user from raw sequence quality control to final phylogenetic inference, making complex bioinformatic pipelines accessible.

## Features

The platform integrates several key bioinformatics tools into a seamless workflow:

-   **Quality Control:** Visualize raw read quality using **FastQC**.
-   **Preprocessing:** Trim adapters and low-quality bases with **Trimmomatic** and merge paired-end reads with **PEAR**.
-   **Sequence Analysis:** Perform OTU (Operational Taxonomic Unit) clustering using **VSEARCH**.
-   **Taxonomic Identification:** Assign taxonomy to OTUs by querying the Barcode of Life Data System (BOLD) using **BOLDigger3**.
-   **Phylogenetic Analysis:**
    -   Perform multiple sequence alignment with **MAFFT**.
    -   Infer a fast, initial phylogenetic tree with **FastTree2**.
    -   Conduct a robust maximum likelihood phylogenetic reconstruction with bootstrap support using **RAxML-NG**.
-   **Interactive Visualization:** View phylogenetic trees directly in the browser.

## Prerequisites & Installation

**This is the most critical step.** The PARAFISH script is cross-platform, but the required external tools must be installed manually according to your operating system.

### 1. Install External Tools

You must install the following software. For Linux/macOS, using a package manager like `apt` or `brew` is highly recommended. For Windows, you will need to download the tools and likely add their containing folders to your system's PATH environment variable.

-   **Java**
    -   **Link:** [java.com](https://www.java.com/en/download/manual.jsp)
    -   **Installation:** Download and install the version for your OS.

-   **Trimmomatic**
    -   **Link:** [usadellab.org](http://www.usadellab.org/cms/index.php?page=trimmomatic)
    -   **Installation:** Download the binary `.zip` file and extract it.

-   **PEAR**
    -   **Link:** [h-its.org](https://www.h-its.org/software/pear-paired-end-read-merger/)
    -   **Installation:** Download the pre-compiled binary for your OS.

-   **FastQC**
    -   **Link:** [bioinformatics.babraham.ac.uk](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    -   **Installation:** Download the source `.zip` file and follow its installation instructions.

-   **VSEARCH**
    -   **Link:** [github.com/torognes/vsearch](https://github.com/torognes/vsearch)
    -   **Installation:** **Linux:** `sudo apt-get install vsearch` | **macOS:** `brew install vsearch` | **Windows:** Download binary from GitHub.

-   **MAFFT**
    -   **Link:** [mafft.cbrc.jp](https://mafft.cbrc.jp/alignment/software/)
    -   **Installation:** **Linux:** `sudo apt-get install mafft` | **macOS:** `brew install mafft` | **Windows:** Download binary.

-   **FastTree**
    -   **Link:** [morgannprice.github.io/fasttree](https://morgannprice.github.io/fasttree/)
    -   **Installation:** **Linux:** `sudo apt-get install fasttree` | **macOS:** `brew install fasttree` | **Windows:** Download binary.

-   **RAxML-NG**
    -   **Link:** [github.com/amkozlov/raxml-ng](https://github.com/amkozlov/raxml-ng)
    -   **Installation:** Follow the compilation or download instructions on their GitHub for your OS.

### 2. Set Up the Python Environment

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/JohnnyFarians24/Project_Bioinformatics.git
    cd Project_Bioinformatics
    ```

2.  **Place Trimmomatic Folder:**
    Ensure the `Trimmomatic-0.39` folder (that you extracted earlier) is placed in the root of the project directory, alongside `PARAFISH.py`.

3.  **Install Python dependencies:**
    It is highly recommended to use a virtual environment.
    ```bash
    python -m venv venv
    # On Windows:
    # venv\Scripts\activate
    # On macOS/Linux:
    # source venv/bin/activate

    pip install -r requirements.txt
    ```

### 3. Configure Tool Paths

This step connects the PARAFISH script to the tools you installed.

1.  In the project folder, you will find a file named `config.ini`.
2.  Open `config.ini` with a text editor.
3.  For each tool listed under `[PATHS]`, provide the full, absolute path to its executable on your system.

**IMPORTANT:**
-   If a tool was successfully added to your system's PATH during installation, you can leave its name as is (e.g., `FASTQC = fastqc`).
-   On **Windows**, always use forward slashes (`/`) or double backslashes (`\\`) for paths.

**Example `config.ini` for a Windows user:**
```ini
[PATHS]
JAVA = C:/Program Files/Java/jdk-11/bin/java.exe
TRIMMOMATIC_JAR = Trimmomatic-0.39/trimmomatic-0.39.jar
PEAR = C:/bio_tools/pear/bin/pear.exe
FASTQC = C:/bio_tools/FastQC/fastqc.bat
...
```

**Example `config.ini` for a macOS/Linux user (if tools are not in PATH):**
```ini
[PATHS]
JAVA = /usr/bin/java
TRIMMOMATIC_JAR = Trimmomatic-0.39/trimmomatic-0.39.jar
PEAR = /Users/myuser/tools/pear/bin/pear
...
```

## Running the Application

Once all prerequisites are installed and the `config.ini` file is correctly set up, run the application from the project's root directory:

```bash
python PARAFISH.py
```

Open your web browser and navigate to `http://127.0.0.1:8050/` to start using the platform.

## Important Notes

### Cross-Platform Compatibility
The PARAFISH script is designed to be compatible with Windows, macOS, and Linux. However, the installation and path configuration of the required external tools is highly specific to your operating system. Please follow the setup instructions carefully.

### Administrator Privileges (`sudo`)
You **do not** need administrator (`sudo`) privileges to **run** the PARAFISH application itself. The application runs with standard user permissions.

You may, however, need administrator privileges to **install** some of the prerequisite tools on Linux and macOS (e.g., `sudo apt-get install ...` or if installing to a system-wide directory).

