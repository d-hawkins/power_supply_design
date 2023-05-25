
6/16/2012

LTSpice has an option to save plots to .wmf however I did not
figure out a method for converting those figures to PDF.

The method I use is;

1) In LTSpice configure the axis to make the figure look nice

2) Configure the Adobe printer in 2-up mode, so the plot only
   takes up half the page (the top half).
   
3) Print to PDF

4) Import into Inkscape
   - LTspice_plot_template.svg has a default setup for
     importing plots, eg., its got a layer with a bounding
     box rectangle, and its got lines where the plot should go.

5) Edit the image
   - Import the image
   - Move it to (10,50) using the X and Y fields at the top
     of the Inkscape GUI
   - Ungroup
   - Delete the page outline and the bottom edge filename
   - Move the title text down by 1mm (so it is closer to the
     correct plot). Check whether the text has a small white
     box behind it, and if it does, delete it.

6) Save as PDF from Inkscape

Plots with lots of lines can take a while to import and manipulate
with Inkscape. The final PDFs look nice though, so its worth being
patient.

LTspice circuits are more of a pain to import. The PDF generated
by Adobe seems to group text in such a way that Inkscape cannot
ungroup it. Rather than trying to edit text in Inkscape, just
delete as much as possible from the LTspice circuit and then
print the 'reduced' version.