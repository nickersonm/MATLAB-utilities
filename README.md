# MATLAB-utilities

Assorted MATLAB utility scripts.


## Description

A set of unrelated MATLAB utility scripts for reuse of common tasks, common functions, and optical simulation.

Some notables:

- [`./Optics/`](./Optics/) contains several free-space scalar beam propagation simulation functions and test scripts.
- [`plotStandard2D`](./plotStandard2D.m) is intended to provide all common and uncommon options desired for plotting 2D lines, including a 2nd y-axis and 2nd x-axis. Read it's documentation to find out more.
- Other `plot*` functions gather common options for various 2D plotting tasks.
- [`figureSize`](./figureSize.m) creates a new figure with a specified handle and size.
- [`fftconv`](./fftconv.m) and [`fftxcorr`](./fftxcorr.m) provide accelerated `conv`olution and `xcorr`elation calculation by transforming to frequency space and multiplying.

> Note: not all functions have been tested with recent MATLAB versions.


## Dependencies

Some functions use the [`smoothn` function](https://www.mathworks.com/matlabcentral/fileexchange/25634-smoothn). This can be replaced with the built-in `smooth` if desired.
