# Truss-Finite-Element-Analysis
A bridge deck is modeled as a combination of bar and beam bending elements. The support is modeled as a standard truss. It is required that we write a code to analyze the structure, and to design the height of the deck and the cross-sectional areas of the truss to make the weight of the bridge as small as possible.

## Run
The code mainly consists of 3 functions.</br>
The first function used is “analyze_structure”. It analyzes the structure given, calculates the stresses in all the members and gets
the ratio between those stresses and the allowable stress.</br>
The second function is “optimize_h” which uses a range of h values and tries different solutions for them using the
“analyze_structure” function. The function then produces the critical h, the ratio, and the stresses in all frame members in that
case.</br>
The third function “optimize_A” uses a range of A values and tries different solutions for them using “analyze_structure”
function. The function then produces the critical A values, the ratio, and the stresses in all truss members in that cases.</br>
The outputs of the run are the height of the frame, the areas of all the truss bars, and the stress in all the structure elements.</br>

### Check the pdf for the problem definition and further details on code implementation.
