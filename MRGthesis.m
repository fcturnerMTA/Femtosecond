clear;
tic

%%%%%%%%%%%%%%%%%%%%
%Defining constants%
%%%%%%%%%%%%%%%%%%%%

%%Assume a value of the Raman gain - figure it out later
g= 1.6e-13;
c= 3e8; %%in m/s

%%The Raman transition
omegaRaman= 775; %%in wavenumbers
dephasingTime= 6.6; %%in picoseconds
HWHMRaman= .2206/(dephasingTime*1e-12)/3e10/2; %%in wavenumbers
%HWHMRaman= 100; %%in wavenumbers
convolution= 0; %%over what TOTAL range should a convolution be carried out in wavenumbers
                 %% Note that a higher value is more accurate, but also much more taxing
                 %% Also note that this will affect the value for g

%%Define the fibre characteristics
FibreLength= 1;  %%in meters
Intervals= 1000;  %%number of points per meter; this number must be large for the program to work
plotPoints= 1000; %%when plotting in 3d, it's best to have only 1000 points
fibreRadius= 250; %%in microns
u11= 2.405; %%zero of the Bessel function J[0](x)
dz= FibreLength/Intervals;

%%Define largest Stokes and Anti-Stokes order you want to include
S_Orders= 20; %%including the long-wavelength beam
A_Orders= 20;

%%You can't interactively change the scale on a surface plot, so set the number you want to plot
minorder= -S_Orders; %%use -S_Orders for all of them
maxorder= A_Orders; %%use A_Orders for all of them

%%Set the frequency resolution
%%Note that a typical bandwidth is about 80 wavenumbers, so a spacing of 16 is appropriate
omegaRes= 7; %%In wavenumbers

%%Set the initial intensity of the Pump, and then the energy ratio between it and the Stokes
Intensity= 6e12; %%in W/cm^2, but I don't think it matters, as g is a parameter to be determined
ratio= 1.5;

%%Finally, set the smallest intensity (relative to the Pump), that you want to plot
plotmin= 1e-4;

%%Labeling convention: short-wavelength beam (Pump beam) is order# S_Orders+3
%%                     long-wavelength beam (Stokes beam) is order# S_Orders+2
%%I'm adding two extra orders on either end for programming considerations
Pump= S_Orders + 3;
Stokes= S_Orders + 2;

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determining resolution and constructing the input spectrum from our data%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omegaIntervals= ceil(omegaRaman/omegaRes); %%number of intervals per order

%%I'm going to make sure that there are an odd number of intervals per order.  By forcing
%%  an odd number, I can define the array elements by using the orders as a reference.
%%  Notice how in the diagram below there are two points on either side of each order
%%  when there's 5 intervals per order.
%%  ..|....|....|....|....|..
%%Equivalently, the # of side elements = (# of intervals - 1)/2

%%The pump frequency is therefore at an index of (Pump Order)*(elements/order)
%%  - the bit on the side

%%  ..|....|....|....|....|....|....|..
%%  1   ^n   ^2n  ^3n  ^Pump*n          the arrows point to the element in question
%%                            n=5, Pump Order=3, S_orders=1, A_orders=1, bit on the side=2


if (mod(omegaIntervals, 2)== 0)
    omegaIntervals= omegaIntervals + 1;
end
side= (omegaIntervals - 1)/2;

%%This now sets the actual resolution, which will be as good as or better than
%%  the specified resolution
omegaRes= omegaRaman/omegaIntervals;

%%Read in the Regen Spectrum
path(path, 'C:/Documents and Settings/Compaq/My Documents/Frasers/Frasers Junk/Waterloo/Work/MRG/Data/Data Analysis/mrg104_308_AND036-AND064');
[lambdaReg, IntensReg]= textread('AND064.mat', '%f %f', 'expchars', 'E', 'commentstyle', 'matlab');
omegaReg= 1e7./lambdaReg;

RegenSpec= zeros(length(omegaReg), 2);
RegenSpec= [omegaReg IntensReg];

%I'm artificially subtracting the background, as I didn't take one to begin with.
RegenSpec(1: length(omegaReg), 2)= sqrt((RegenSpec(1: length(omegaReg), 2) - 9e-7).^2);
RegenSpec(1: length(omegaReg), 2)= sqrt((RegenSpec(1: length(omegaReg), 2) - 3.5e-7).^2);
RegenSpec(1: length(omegaReg), 2)= sqrt((RegenSpec(1: length(omegaReg), 2) - 1.75e-7).^2);
RegenSpec(1: length(omegaReg), 2)= sqrt((RegenSpec(1: length(omegaReg), 2) - 0.9e-7).^2);
RegenSpec(705: length(omegaReg), 2)= 0;
RegenSpec(220: 620, 2)= 0;
RegenSpec(1: 70, 2)= 0;

%%Determine the frequency spacing
points= length(omegaReg);
range= omegaReg(1) - omegaReg(points);
dOmegaReg= range/points;

%%Determine the desired index spacing from the specified resolution.  This should
%%  give a resolution that is as good or better than that specified
spacing= floor(omegaRes/dOmegaReg);
omegaRes= spacing*dOmegaReg;

RegIntervals= floor(points/spacing); %%number of intervals over entire Regen range

%%This is an example of what the array could look like; this time the lines are
%%  showing the array elements we want to keep
%%    |........|........|........|....  "points= 32, spacing= 9, intervals= 3"
%%   1^     9+1^   2*9+1^   3*9+1^ 32^      Listing the index at the arrow

%%Now to get rid of all the extraneous points.  Follow the diagrams above and below
%%    |........|........|........|....
%%    |........|........|........|
%%    ||........|........|
%%    |||........|
%%    ||||
RegenSpec((spacing*RegIntervals + 2): points, :)= [];
for n= 1: RegIntervals
    RegenSpec(1 + n: spacing + n - 1, :)= [];
    omegaReg(1 + n: spacing + n - 1, :)= [];
end

%%Now to set the ratio between the pump and the Stokes.  We measure energy, so integrate
%%  the two peaks in the spectrum and determine the ratio of the areas.
[pumpmax, pumpindex]= max(RegenSpec(:, 2));
omegaPump= omegaReg(pumpindex);
RegenSpec(:, 2)= RegenSpec(:, 2)./pumpmax; %normalize the array to the Pump intensity

StokesArea= 0;
PumpArea= 0;

%%In the original file, the Stokes beam is between the 70th and 220th data point,
%%  and the Pump beam is between the 570th and 710th data point   
for n= floor(70/spacing): ceil(220/spacing)
    StokesArea= StokesArea + (RegenSpec(n, 2) + RegenSpec(n - 1, 2))/2;
end

for n= floor(570/spacing): ceil(710/spacing)
    PumpArea= PumpArea + (RegenSpec(n, 2) + RegenSpec(n - 1, 2))/2;
end

SpectraRatio= PumpArea/StokesArea;
scale= ratio/SpectraRatio;

RegenSpec(floor(70/spacing): ceil(220/spacing), 2)= RegenSpec(floor(70/spacing): ceil(220/spacing), 2)...
                                                        ./scale;

%%Initialize the Frequency and Intensity arrays
%%I'm adding a couple extra orders on either end as it's needed later
omega= (omegaPump - (omegaIntervals*(S_Orders + 2) + side)*omegaRes: omegaRes: ...
    omegaPump + (omegaIntervals*(A_Orders + 2) + side)*omegaRes);
%%Translation: Start at the centre, go back by ((the number of Stokes orders + two extra)
%%  *(number of array elements per order) plus the bit on the side),
%%  increment by the frequency resolution until you get to (Pump + Anti-Stokes +
%%  two extra orders), then add the number of elements on the side

In= zeros(1, length(omega));
for n= pumpindex: length(RegenSpec(:, 2));
    In(Pump*omegaIntervals - side - n + pumpindex)= RegenSpec(n, 2);
end

for n= 1: pumpindex - 1;
    In(Pump*omegaIntervals - side + n)= RegenSpec(pumpindex - n, 2);
end

%GVD= 0000e-6; %%I've set it to 10000; for Steve it was -1877
%TOD= 000000e-9; %%I've set it to 1000000; for Steve it was 18719
%extraphase= GVD/2.*(omega*100*3e8*2*pi).^2 + TOD/6.*(omega*100*3e8*2*pi).^3;

%lambda= 1e4./omega; %%this gives the wavelength in microns

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the phase dispersion%
%%%%%%%%%%%%%%%%%%%%%%%%%%

ngas= zeros(1, length(omega));
k= zeros(1, length(omega));
%ngas= ngas + 1; %%comment what's below to set no dispersion

%%Sellmeier equations for the gas
%%THIS IS FOR CO2 AT STP, but I need to use something so I'll start with this
%ngas= 1 + 0.012*(0.58*lambda.^2./(166*lambda.^2 - 1) + 0.12*lambda.^2./(79.6*lambda.^2 - 1)...
%    +0.0053*lambda.^2./(56.3*lambda.^2 - 1) + 0.00432*lambda.^2./(46*lambda.^2 - 1)...
%    +0.000122*lambda.^2./(0.058*lambda.^2 - 1));

ngas= 1 + 0.7e30./((31.9e15)^2 - (2.*pi.*c.*100.*omega).^2); %%from physrevA68_023812p15_2003

%%Sellmeier equation for the silica fibre
%nsilica= sqrt(1.2955 + 0.80985*lambda.^2./(lambda.^2 - .0107945) + 0.91714*lambda.^2./(lambda.^2 - 100));
nsilica= 1.5;
n1= (nsilica.^2 + 1)./(2.*sqrt(nsilica.^2 - 1));

%%Axial wavevector for the EH11 hybrid mode
k= 2.*pi.*omega.*100.*ngas.*(1 - 1/2.*(u11./(omega.*100)./(2.*pi.*fibreRadius)).^2).*(1 - i*n1./(omega.*100)/pi/fibreRadius);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Electric Field Amplitudes%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc;
A= zeros(size(In));
B= zeros(size(In));
C= zeros(size(In));
P= zeros(1, Intervals + 1);

A= sqrt(In); 
B= A;
C= A;
%%I'll use A, B, and C to keep track of things as I go through the program

%%These next two things are what I'll use to graph the pulses
%surfacegraph= zeros(Intervals+1, (A_Orders + S_Orders + 5)*omegaIntervals);

%%Propagate in z
for z= 0: Intervals
    
    %%Sum over the absorption(i)->emission(j) events within a reasonable range
    %%  Convolving each contribution with a particular detuning range
    for x= -ceil(convolution/2/omegaRes): ceil(convolution/2/omegaRes)
        %%I don't want to include the extra orders I put in, so start just beyond the
        %%  first order and go to the end of the Anti-Stokes (see diagram below)
        
        %%  ..|....|....|....|....|....|....|..
        %%  1   ^n   ^2n  ^3n  ^Pump*n          the arrows point to the element in question
        %%                            n=5, Pump Order=4, S_orders=1, A_orders=1, bit on the side=2

        %%Initialize the material excitation for each detuning
        %%  (de-excitation is given by its complex-conjugate)
        Q= 0;
                
        for n= 3*omegaIntervals + 1: (Pump + A_Orders)*omegaIntervals
        %%Using the definition of Q below, I need to start from one order above the mininmum
        
%            Omega= omega(n) - omega(n - omegaIntervals + x);
            Q= Q + B(n)*conj(B(n - omegaIntervals + x))...
                *exp( -i*(k(n) - k(n - omegaIntervals + x))*dz);%...
%                *1/(omegaRaman^2 - Omega^2 - 2*i*Omega*HWHMRaman);
        end
        
        P(z + 1)= real(Q);
        
        %%Increment each order
        for n= 2*omegaIntervals + 1: (Pump + A_Orders)*omegaIntervals
            C(n)= B(n) + g*Intensity/ngas(n)*omega(n)/omega(Pump)*...
                (conj(Q)*B(n + omegaIntervals - x)*exp(i*(k(n) - k(n + omegaIntervals - x))*dz)...
                     - Q*B(n - omegaIntervals + x)*exp(i*(k(n) - k(n - omegaIntervals + x))*dz))*dz;
        end
    end    

    %%Now to define the values for a 3d plot in frequency and distance
    for n= 2*omegaIntervals + 1: (Pump + A_Orders)*omegaIntervals
        surfacegraph(z + 1, n)= C(n).*conj(C(n)).*Intensity;
    end
    
    B= C;
end


%%%%%%%
%Plots%
%%%%%%%

toc;
%%Now I'm changing the convention - the Pump is of order 0, and the Stokes is of order -1
X= ( -(S_Orders + 2) - side/omegaIntervals: 1/omegaIntervals: (A_Orders + 2) + side/omegaIntervals);
Y= ((C.*conj(C))*Intensity);
L= log10((C.*conj(C) + plotmin)*Intensity);
figure(1)
%plot(X, Y)
plot(X, L)
xlim([-0.5 7])

SX= (0: FibreLength/Intervals: FibreLength);
SY= (minorder - side/omegaIntervals: 1/omegaIntervals: maxorder + side/omegaIntervals); %%See diagram below

%%    ..|....|....|....|....|....|....|..
%%               min   P   max

%%The 3d plot works okay for 500 intervals; above and below that it gets messy
%plotinterval= floor(Intervals/500);

    %%Delete elements up to just before the side of the minorder (use the diagram)
    surfacegraph(:, 1: (Pump + minorder - 1)*omegaIntervals - 1) = [];
    
    %%Delete elements beyond the side of the maxorder, note that I've already deleted
    %%  some elements in the line above (the 2 is because of that extra order I added earlier)
    surfacegraph(:, (-minorder + maxorder + 1)*omegaIntervals:...
        (-minorder + A_Orders + 3)*omegaIntervals - 1)= [];

    %%Delete elements in z to have only 1000 elements - this makes the best plot
    plotInterval= floor((Intervals)/plotPoints);
    for n = 1: plotPoints - 1
        surfacegraph((n + 1): (n + plotInterval - 1), :) = [];
        SX((n + 1): (n + plotInterval - 1)) = [];
    end
    if (length(surfacegraph(:, 1)) > plotPoints)
        surfacegraph((plotPoints + 1): length(surfacegraph(:, 1)), :) = [];
    end
    if (length(SX) > plotPoints)
        SX((plotPoints + 1): length(SX)) = [];
    end
SZ= surfacegraph;
    %%Delete some elements so that there's at most 500 points along z
%    if (Intervals> 500)
%        for n= Intervals + 1: -1: 1
%            if (mod(n, plotinterval)~= 1)
%                SZ(n, :)= [];
%            end
%        end
%    end

%SL= log10(plotmin.*Intensity + surfacegraph); 
    
%SL = loggraph;
%    SL(:, 1: (Pump + minorder - 1)*omegaIntervals - 1) = [];
%    SL(:, (-minorder + maxorder + 1)*omegaIntervals:...
%        (-minorder + A_Orders + 3)*omegaIntervals) = [];
%    if (Intervals> 500)
%        for n= Intervals + 1: -1: 1
%            if (mod(n, plotinterval)~= 1)
%                SL(n, :)= [];
%                SX(n)= [];
%            end
%        end
%    end
    
figure(2)
mesh(SY, SX, SZ)
%mesh(SY, SX, SL)
%shading interp
toc