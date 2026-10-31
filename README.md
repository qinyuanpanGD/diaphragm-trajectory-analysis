# Diaphragm trajectory analysis

Consensus k-means using DE, DTee and DTF on days 1, 4 and 7. Defaults are 1,000 resampling iterations, 80% participant sampling and K = 2 to 6.

Input CSV columns: `ID`, `time`, `DE`, `DTee`, `DTF`.

```sh
Rscript install_packages.R
Rscript analysis.R "data/analysis.csv" "results"
```

Participant-level study data are not included. Output files are saved in the specified results directory.
