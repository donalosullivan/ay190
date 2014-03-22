<h2>Project Description</h2>

<p>
This <i>Ay190: Computational Astrophysics'</i> project takes visibilities from a radio interferometry array and constructs an image of the source using both a Discrete Fourier Transform method and NumPy's Fast Fourier Transform method.
</p>

<h2>Input Files</h2>

<i>'AntennaPositions.csv'</i> - This is a text file that contains the position (of each antenna) measured in centimeters. This first column gives the antenna number, the second column gives the x-position of the antenna, and the final column gives the y-position of the antenna.

<i>'Visibilities.csv'</i> - This file gives the visibilities measured by the radio interferometer. These visibilities have no thermal noise and are perfectly phase-calibrated (but they are not flux-calibrated so the amplitude is meaningless). The first two columns define the antennas that the visibility was measured between (for example, 2 and 17 indicates that the visibility was measured between antennas 2 and 17). The final two columns give the amplitude and phase of the visibility respectively (the complex visibility is then amplitude*exp(i*phase) ).

These observations have been made at a wavelength of 1cm (with an infinitesimally small bandpass filter so that there is no bandwidth smearing). Additionally, the source is located directly overhead at the zenith (so there are no projection effects), and the beam of each antenna is Gaussian with sigma=1 arcminute. The integration time is too small for the rotation of the Earth to have any impact on the source's position in the sky.


<h2>Data Structures</h2>

<h3><b><i>'pos'</i> - Antenna Positions</b></h3>
    pos( (i,x,y) )<br /> 
        i is the antenna number, x and y are cartesian coordinates<br /> 
        -->i,x,y = pos(0) would unpack the number and position of the first antenna<br /> 

<h3><b><i>'vis'</i> - Positions and Visibilities</b></h3>
    vis( (i,j,A,phi) )<br /> 
        i and j are indices of the two antennae in 'pos', <br /> 
        A is amplitude,<br /> 
        phi is phase <br /> 
    
<h3><b><i>'uvvis'</i> - Baselines and Visibilities </b></h3>
    uvvis( (u,v,A,phi) )<br /> 
        for each row in 'vis', 'uvvis' stores:<br /> 
            u,v - the baseline vectors (x2-x1) and (y2-y1) of antennae i and j<br /> 
            A,phi - the same visibility as in 'vis'<br /> 

<h3><b><i>'l' and 'm'</i> - RA and DEC grids</b> </h3>
    These 1D arrays are grids of RA and DEC ranging from -100'' to +100''<br /> 
    Centered on zenith<br /> 



<h2>Methods</h2>

<h3><b>A(l,m,sig=arcmin)</b></h3>
Gaussian function for Antenna Beam A(l,m)<br /> 

<h3><b>DFT_rhs(uvvis,_l,_m)</b></h3>
RHS function for Discrete FT: Sum over (u,v)<br />
   
    
<h3><b>DFT(uvvis,L,M)</b></h3>
Returns intensity array I(l,m) using DFT_rhs<br />
    
    
<h3><b>find_nearest_gridpoint(x,xlist)</b></h3>
Takes a measured (u,v) point and locates the nearest (u,v) point on an evenly spaced grid<br />
    

<h3><b>uv_grid(uvvis,ugrid,vgrid)</b></h3>
Create an evenly spaced grid of visibilities to use in an inverse fft<br />
  

<h3><b>def get_selection(pos,vis,N,orderby='asc')</b></h3>
Return indices of rows in uvvis which only use N closest/farthest antennae<br />






