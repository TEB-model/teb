# TEB Tests

This folder contains regression and integrations tests used to develop TEB. All tests are based on data from the CAPITOUL experiment.

All tests are executed using with the following syntax


path_to_exe_ref = build_teb(commit_id_ref)
path_to_exe_this = build_teb(commit_id_this)
compare_teb(path_to_exe_ref, path_to_exe_this, threshold=0, plots=True)