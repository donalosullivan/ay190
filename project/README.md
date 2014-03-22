<h2>Project Description</h2>

<p>
This <i>Ay190: Computational Astrophysics'</i> project takes visibilities from a radio interferometry array and constructs an image of the source using both a Discrete Fourier Transform method and NumPy's Fast Fourier Transform method.
</p>

<h2>Input Files</h2>

<i>'AntennaPositions.csv'</i> - This is a text file that contains the position (of each antenna) measured in centimeters. This first column gives the antenna number, the second column gives the x-position of the antenna, and the final column gives the y-position of the antenna.

<i>'Visibilities.csv'</i> - This file gives the visibilities measured by the radio interferometer. These visibilities have no thermal noise and are perfectly phase-calibrated (but they are not flux-calibrated so the amplitude is meaningless). The first two columns define the antennas that the visibility was measured between (for example, 2 and 17 indicates that the visibility was measured between antennas 2 and 17). The final two columns give the amplitude and phase of the visibility respectively (the complex visibility is then amplitude*exp(i*phase) ).

These observations have been made at a wavelength of 1cm (with an infinitesimally small bandpass filter so that there is no bandwidth smearing). Additionally, the source is located directly overhead at the zenith (so there are no projection effects), and the beam of each antenna is Gaussian with sigma=1 arcminute. The integration time is too small for the rotation of the Earth to have any impact on the source's position in the sky.


<h2>Data Structures</h2>

<i>'pos'</i> - Antenna Positions <br />
    pos( (i,x,y) )<br /> 
        i is the antenna number, x and y are cartesian coordinates<br /> 
        -->i,x,y = pos(0) would unpack the number and position of the first antenna<br /> 

<i>'vis'</i> - Positions and Visibilities <br />
    vis( (i,j,A,phi) )<br /> 
        i and j are indices of the two antennae in 'pos', <br /> 
        A is amplitude,<br /> 
        phi is phase <br /> 
    
<i>'uvvis'</i> - Baselines and Visibilities <br />
    uvvis( (u,v,A,phi) )<br /> 
        for each row in 'vis', 'uvvis' stores:<br /> 
            u,v - the baseline vectors (x2-x1) and (y2-y1) of antennae i and j<br /> 
            A,phi - the same visibility as in 'vis'<br /> 

<i>'l' and 'm'</i> - RA and DEC grids <br />
    These 1D arrays are grids of RA and DEC ranging from -100'' to +100''<br /> 
    Centered on zenith<br /> 



<h2>Methods</h2>

<b>A(l,m,sig=arcmin)</b><br />
Gaussian function for Antenna Beam A(l,m)<br /> 

<b>DFT_rhs(uvvis,_l,_m)</b><br />
RHS function for Discrete FT: Sum over (u,v)<br />
   
    
<b>DFT(uvvis,L,M)</b><br />
Returns intensity array I(l,m) using DFT_rhs<br />
    
    
<b>find_nearest_gridpoint(x,xlist)</b><br />
Takes a measured (u,v) point and locates the nearest (u,v) point on an evenly spaced grid<br />
    

<b>uv_grid(uvvis,ugrid,vgrid)</b><br />
Create an evenly spaced grid of visibilities to use in an inverse fft<br />
  

<b>def get_selection(pos,vis,N,orderby='asc')</b><br />
Return indices of rows in uvvis which only use N closest/farthest antennae<br />






