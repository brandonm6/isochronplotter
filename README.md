# IsochronPlotter
Constructs isochron diagrams for argon-argon (40Ar/39Ar) data.

## Method
1. Iterates down from 1/298.6 (the 36Ar/40Ar ratio of atmospheric argon) using increments of 99% of the previous ratio 
   1. stops at the lowest 36Ar/40Ar ratio in the split of the sample 
2. Recalculates the age of each step using a 40Ar*/39Ar, calculated with the current 36Ar/40Ar ratio 
3. Finds groups of at least three consecutive steps (consecutive steps will have an age interval that overlaps) and 
checks if they have at least 50% cumulative 39Ar - these are considered plateaus 
4. Draws trendlines using the groups of steps identified as plateaus (uses regression based on the approach by Mahon 1996)

## Usage
### Preparing a data set
IsochronPlotter accepts CSVs formatted with the following columns:

| Run_ID | Isoch_39_Over_40 | Pct_i39_Over_40_Er | Isoch_36_Over_40 | Pct_i36_Over_40_Er | Correl_36_Over_39 | Ar36_Over_Ar39 | Ar39_Moles | Age_Er | J |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |---|
<br>
Full raw run data downloaded from Mass Spec in CSV format should already contain these columns.

### Plotting an isochron diagram
1. Place your appropriately formatted CSV file in the _csvs_ folder.
2. In _run.py_, change `file_name` to the name of your CSV.
3. Run _run.py_.

ex.
```
file_name = "CIA dense 197'.csv"
```

PNGs of the plots are saved to `isochronplotter/outputs`.

#### Attributes of Plotter

- omit: list of letters or numbers to remove from all splits (ex. `["a", "232226-01", "1"]`)
- save: download plots as PNGs to `isochronplotter/outputs`
- verbose: print initial <sup>40</sup>Ar/<sup>36</sup>Ar ratios that have at least three consecutive steps and all 
plateaus found by LocatePlateaus (_locateplateaus.py_)
- iverbose: print DataFrames for each increment along with the groups of overlapping steps

#### Useful variables to call from Plotter

- `force_remove`: list of removed run-ids that had NaN values
- `df`: DataFrame holding the modified data

### Prerequisites
1. Numpy
2. Pandas
3. Matplotlib

---
For questions, contact Brandon Maher at brandmaher@gmail.com
