# MC generator comparisons

This readme discribes how to run the post analysis on the the `AnalysisResults.root` file produced by `/PWG/HMTF/macros/AddTaskHMTFMCMultEst.C`

## Dependencies

The post analysis is implemented in python and depends on the `rootpy` (www.rootpy.org) module and `LaTex`. The easiest and cleanest way to maintain a python develpment environment is by installing `virtualenvwrapper` and using the python package installer `pip`. `virtualenvwrapper` is like Dario's aliroot script for python modules

	$ sudo apt-get install python-pip
	$ pip install --user virtualenvwrapper # this installs virtualenvwrapper to ~/.local
	$ mkvirtualenv mcstudies   # Create the environment `mcstudies`
	$ workon mcstudies    # source the environment `mcstudies; use '$ deactivate' to leave the virtual environment
	$ pip install rootpy  # install rootpy into the mcstudies environment

## Running the post analysis

The post analysis produces a summary pdf and also write all the plots to the given input file to the folder `Results{trigger}`. To run the analysis do:

	$ workon mcstudies   # if its a new shell
	$ python ./post_main.py path/to/AnalysisResults.root {Inel, InelGt0, V0AND} "<summary name>"

Where Inel, InelGt0 or V0AND selects the trigger for which the Summary should be produced.

This produces a new folder in the current working directory derived from the summary title. Within this folder are is the summary pdf file.
