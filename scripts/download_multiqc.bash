# Download multiqc data
DATADIR=/gscratch/scrubbed/elpetrou/pollock/fastqc/Lane2
MYDATA1=multiqc_report.html
MYDATA2=multiqc_data
TARGETDIR=C:\Users\elpet\OneDrive\Documents\pollock

# Download multiqc html file
scp elpetrou@klone.hyak.uw.edu:$DATADIR'/'$MYDATA1 $TARGETDIR

# download folder with plot data for multiqc
scp -r elpetrou@klone.hyak.uw.edu:$DATADIR'/'$MYDATA2 $TARGETDIR
scp -r elpetrou@klone.hyak.uw.edu:/gscratch/scrubbed/elpetrou/pollock/fastqc/Lane2 ~

# using putty
# In start, type cmd to get a windows command window
set PATH="%PATH%;%ProgramFiles%\putty"

pscp elpetrou@klone.hyak.uw.edu:/gscratch/scrubbed/elpetrou/pollock/fastqc/Plate3/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt C:\Users\elpet\OneDrive\Documents\pollock