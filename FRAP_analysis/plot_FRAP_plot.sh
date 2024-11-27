#example of how to run plotting script, with input of xlsx file with the following columns:
#"FRAP" (a unique identifier for each cell/site), "Time", "condition", "BG" (background intensity measurements to subtract), "Mean" (mean intensity for the ROI)

#Fig 1: DYRK3, GFP, PPP1CC
python plot_frap_results.py --xlsx DYRK3_PP1_combined_frap.xlsx -s SRRM2-mCh SRRM2-mCh+GFP SRRM2-mCh+GFP-DYRK3 SRRM2-mCh+YFP-PPP1CC -c \#000000 \#2E3532 \#ae7100 \#2626AC -n fig1_frap
