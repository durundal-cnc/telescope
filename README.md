This program drives an alt/az telescope mount to a variety of targets, including astronomical features, and satellites tracked by the NORAD catalogue. It uses an iPhone to record the current location and orientation as a starting point, correction from observation of astronomical features is not yet implemented. Currently an internet connection is required to get the two line ephemeris (TLE) from spacetrack.org, though the actual trajectory computations are all offline and requisite TLEs could be downloaded ahead of time if going somewhere without signal.

The hardware I'm using is driven with a pair of harmonic drive servomotors for each axis, coupled to a reduction timing belt to yield appropriate angular resolution. You can adapt the gear reduction parameters for your specific system, and add backlash compensation if needed (currently set to 0 thanks to the harmonic drives).

It also includes support for taking a grid of photos to make a panoramic stitched photograph, by specifying starting and ending coordinates and then mapping the space between them.

In work as a fun hobby project! GUI layout still very much in work.