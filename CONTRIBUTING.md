# Contributing to svgpathtools

The following is a few and guidelines regarding the current philosophy, style, 
flaws, and the future directions of svgpathtools.  These guidelines are meant 
to make it easy to contribute.

## Being a Hero
We need better automated testing coverage.  Please, submit unittests!  See the 
Testing Style section below for info.

Here's a list of things that need (more) unittests:
* TBA (feel free to help)

## Submitting Bugs
If you find a bug, please submit an issue along with an **easily reproducible 
example**.  Feel free to make a pull-request too (see relevant section below).


## Submitting Pull-Requests 

#### New features come with unittests and docstrings.
If you want to add a cool/useful feature to svgpathtools, that's great!  Just 
make sure your pull-request includes both thorough unittests and well-written 
docstrings.  See relevant sections below on "Testing Style" and 
"Docstring Style" below.


#### Modifications to old code may require additional unittests.
Certain submodules of svgpathtools are poorly covered by the current set of 
unittests.  That said, most functionality in svgpathtools has been tested quite 
a bit through use.
The point being, if you're working on functionality not currently covered by 
unittests (and your changes replace more than a few lines), then please include 
unittests designed to verify that any affected functionary still works.


## Style 

### Coding Style
* Follow the PEP8 guidelines unless you have good reason to violate them (e.g. 
you want your code's variable names to match some official documentation, or 
PEP8 guidelines contradict those present in this document).
* Include docstrings and in-line comments where appropriate.  See 
"Docstring Style" section below for more info.
* Use explicit, uncontracted names (e.g. `parse_transform` instead of 
`parse_trafo`).   Maybe the most important feature for a name is how easy it is 
for a user to guess (after having seen other names used in `svgpathtools`).
* Use a capital 'T' denote a Path object's parameter, use a lower case 't' to 
denote a Path segment's parameter.  See the methods `Path.t2T` and `Path.T2t` 
if you're unsure what I mean.  In the ambiguous case, use either 't' or another 
appropriate option (e.g. "tau"). 


### Testing Style
You want to submit unittests?!  Yes!  Please see the svgpathtools/test folder 
for examples.


### Docstring Style
All docstrings in svgpathtools should (roughly) adhere to the Google Python 
Style Guide.  Currently, this is not the case... but for the sake of 
consistency, Google Style is the officially preferred docstring style of 
svgpathtools.  
[Some nice examples of Google Python Style docstrings](
https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html)
