RA3EventDisplay
===============

Simple event display for RA3 ntuples using the PF particle list

Functions
Basic function is to display events in eta-phi plane. Each object is colored according to
the particle flavor assigned by the PF algorithm, and has the size proportional to the log
of its Pt. One can choose to work with different PFJet collections (default sets stored in
the ntuples being ak5 and ak5chs), or with multiple MET collections.
There are interfaces to calculate the invariant mass and the Mt of a given set of particles,
but this is not quite polished in terms of usability (see below).

Basic usage
 RA3EventDisplay evdisp;
 evdisp.addPath("path_to_susyEvents.root"); // can use wildcard to add multiple files
 evdisp.showEvent(run, event); // seeks to the given event number
 evdisp.showNextEvent(); // increments entry number by one

IMPORTANT NOTE
RA3 ntuples store a jet only when its L1FastL2L3-corrected Pt is above a threshold (default
20 GeV). Similarly, PF particle is stored only if its Pt is above a threshold (default 3 GeV),
or if it is a part of at least one jet collection.
This setup might result in strange appearance of the clusters of objects in the display. For
example, if a cluster belongs to an ak5chs jet but doesn't to an ak5 jet, then when using the
ak5 for the display, you will see non-jet clusters consisting of very soft particles. It is
important to understand that what is displayed are a somewhat arbitrarily chosen part of the
whole event.

Changing the parameters of the display
setPtThreshold(double pt): sets the lower limit on the Pt of the particles to be displayed
setPFJetCollection(TString const& col): sets the name of the jet collection to be displayed
addMet(TString const& met, bool add = true, int color = 0): adds a MET named met if add is
 true. If false, removes the corresponding met. If color != 0, a dashed line of the given
 color is used to represent the MET.
setMainMet(TString const& met): when multiple METs are displayed, sets the main MET whose
 magnitude is printed at the top right corner of the display. Main MET is also used to
 calculate the Mt (see below).


Invariant mass and Mt calculations
Each PF particle in the display is represented by a TMarker, which cannot be named. As such,
the identification of the particles are done by using the fUniqueID field of the TObject. To
find the unique ID, you need to right-click the marker of interest and select "Inspect" from
the pull-down menu.
The invariant mass and Mt calculation functions take this unique ID (an integer) as the
argument. For Mt calculation, the "main MET" is used.
