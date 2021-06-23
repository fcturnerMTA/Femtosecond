clear all;
close all;
tic

%%%%%%%%%%%%%%%%%%%%
%Defining constants%
%%%%%%%%%%%%%%%%%%%%

halt= 0;
%%Assume a value of the Raman gain - figure it out later
g= 6e-12;

%%Define the interaction length and interval spacing
Fibre_Length= 1;  %%in meters
fibreRadius= 64e-6; %%in meters
c= 3e8; %%in m/s
u11= 2.405; %%hollow-fibre mode constant (zero of J0(u)) for the EH11 mode
Pressure= 1.5; %%in atm

Intervals= 3000000;
plotPoints= 1000;  %%number of slices along the fibre, 100 is typical
dz= Fibre_Length/Intervals;

omegaPump= 12500; %%in wavenumbers
omegaRaman= 775; %%in wavenumbers

%%Define largest Stokes and Anti-Stokes order you want to include
S_Orders= ceil(omegaPump/omegaRaman) - 1; %%including the long-wavelength beam
A_Orders= 200;

%%Labeling convention: short-wavelength beam (Pump beam) is order# S_Orders+2
%%                      long-wavelength beam (Stokes beam) is order# S_Orders+1
%%I'll be adding an extra order on each end, which is why I add the extra 1
Pump= S_Orders + 2;
Stokes= S_Orders + 1;

%%Set the initial intensity of the Pump, and then the ratio between it and the Stokes
Intensity= 6e12; %%in W/cm^2, but I don't think it matters, as g is a parameter to be determined
ratio= 1;

%%Finally, set the smallest intensity (relative to the Pump), that you want to plot
plotmin= 10^(-4.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Electric Field Amplitudes%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Set the maximum frequency spacing
%%For now I'm setting this to simply resolve individual plane-wave peaks
omegaRes= omegaRaman/3; %%In wavenumbers

%%Now to figure out how large each order should be. I'm rounding up, which
%%  will give a resolution as good or better than that specified above
number= ceil(omegaRaman/omegaRes);
if (mod(number, 2)== 0)
    number= number + 1;
end
omegaRes= omegaRaman/number;

%%                      ..|....|....|....|....|..
%%By setting an odd number of spaces, I can define the number of elements in the
%%  arrays by using the orders as a reference.  Notice how in the above picture
%%  there are two points on either side of each line when there's 5 spaces.
%%Equivalently there will be (number - 1)/2 array elements on either side of each order.
side= (number - 1)/2;

%%Initialize the Frequency, Wavevector, Intensity, and EField Amplitude arrays
%%I'm adding an extra order on either end as it's needed later
omega= (omegaPump - (number*(S_Orders + 1) + side)*omegaRes: omegaRes: ...
    omegaPump + (number*(A_Orders + 1) + side)*omegaRes);
%%Translation: Start at the centre, go back by the (number of Stokes orders + one extra)
%%  *(number of array elements per order) plus the bit on the side),
%%  increment by the frequency resolution until you get to the Pump + Anti-Stokes +
%%  one extra order, then add the number of elements on the side
nglass= 1.5;
n1= 0.5*(nglass^2 + 1)/sqrt(nglass^2 - 1);
ngas= sqrt(1 + Pressure.*0.7e30./((31.9e15)^2 - (2.*pi.*c.*100.*omega).^2)); %%from physrevA68_023812p15_2003
%ngas= 1;
lambda= 1./(omega.*100.*ngas);
k= 2.*pi./lambda.*(1 - 1/2.*(u11./(omega.*100)./(2.*pi.*fibreRadius)).^2.*(1 - i.*n1./(100.*omega.*pi.*fibreRadius)));
%alpha= (u11./2./pi.*lambda).^2./fibreRadius.^3.*n1;
%alpha= 0;

%%Setting the intensity of the cw pump and stokes
In= zeros(1, number*(S_Orders + A_Orders + 2) + 2*side + 1);
In(Pump*number - side)= Intensity;
In(Stokes*number - side)= Intensity*ratio;
%%IMPORTANT!!!!  The pump frequency is now at an index of (Pump Order)*(elements/order)
%%  - the bit on the side

%%  ..|....|....|....|....|..
%%  1   ^n   ^2n  ^Pump*n          the arrows point to the element in question
%%                            n=5, Pump Order=3, S_orders=1, A_orders=1, bit on the side=2

A= sqrt(In)/sqrt(Intensity);
B= A;
C= A;
%%I'll use A, B, and C to keep track of things as I go through the program

%%These next two things are what I'll use to graph the pulses
%surfacegraph= zeros(Intervals+1, (A_Orders + S_Orders + 3)*number);
%loggraph= (zeros(Intervals+1, (A_Orders + S_Orders + 3)*number)) + log10(plotmin*Intensity);
%%Using the plotmin*Intensity to avoid log(zero)

toc;
%%Propagate in z
for z= 0: Intervals
    %%Initialize the material excitation (de-excitation is given by its complex-conjugate)
    Q= 0;
    
    %%Sum over the absorption(i)->emission(j) events that give omega(i)-omega(j)= Omega
    %%Off-resonant contributions won't come into effect without larger bandwidths
    for n= (Pump - S_Orders)*number: (Pump + A_Orders)*number
        %%I don't want to include the extra orders I put in, so start just beyond the
        %%  first order and go to the end of the Anti-Stokes (see diagram a ways up)
        Q= Q + Intensity*B(n)*conj(B(n - number));%*exp( -i*(k(n) - k(n - number))*dz);
    end
    
    %%Keep track of Q with z
    %P(z + 1)= real(Q);
    
    %%Increment each order
    for n= number + 1: (Pump + A_Orders)*number
        C(n)= B(n) + g*Pressure*omega(n)/omega(Pump*number - side)*...
            (conj(Q)*B(n + number)*exp(i*(k(n) - k(n + number))*dz)...
              - Q*B(n - number)*exp(i*(k(n) - k(n - number))*dz))*dz;
        
        C(n)= C(n).*exp(i.*k(n).*dz);
    end

%    if (halt== 0)
%        D(:, z + 1)= C(:);
%        F(:, z + 1)= fft(D(:, z + 1));
%        F(:, z + 1)= fftshift(F(:, z + 1));
                
%        if (z~= 0)
            if (mod(z, floor(Intervals/plotPoints))== 0)
                D(:, z*plotPoints/Intervals + 1)= C(:);
                F(:, z*plotPoints/Intervals + 1)= fft(D(:, z*plotPoints/Intervals + 1));
                F(:, z*plotPoints/Intervals + 1)= fftshift(F(:, z*plotPoints/Intervals + 1));
%                halt= 1;
            end
%        end
%    end
    B= C;
end

toc

%%%%%%%
%Plots%
%%%%%%%

%%Now I'm changing the convention - the Pump is of order 0, and the Stokes is of order -1
X= ( -(S_Orders + 1) - side/number: 1/number: (A_Orders + 1) + side/number);
Y= ((C.*conj(C))*Intensity);
%plot(X, Y)

figure(1)
plot(X, log10(Y + plotmin*Intensity), 'LineWidth', 2)
xlim([-18 25])

figure(2)
realD= real(D)./max(max(real(D)));
InFrequency= D.*conj(D)./(max(max(D)).*conj(max(max(D))));

pause
aviobj= avifile('spec.avi', 'fps', 60);
for n= 1: plotPoints + 1
    h= plot(X, log10(InFrequency(:, n) + plotmin), 'LineWidth', 2);
    xlim([-20 40])
    ylim([-4.5 0])
    set(h,'EraseMode','xor');
    frame= getframe(gca);
    aviobj= addframe(aviobj,frame);
end
aviobj= close(aviobj);

figure(3)
realF= real(F)./max(max(real(F)));
InTime= F.*conj(F)./(max(max(F)).*conj(max(max(F))));
time= (1:length(C))./omegaRes./length(C)./3e10.*1e15;
time= time - time(round(length(time)/2));

aviobj= avifile('fft.avi', 'fps', 60);
for n= 1: plotPoints + 1
    h= plot(time, InTime(:, n), 'LineWidth', 2);
    ylim([0 1]);
    set(h, 'EraseMode', 'xor');
    frame= getframe(gca);
    aviobj= addframe(aviobj,frame);
end
aviobj= close(aviobj);

figure(4)
aviobj= avifile('efield.avi', 'fps', 60);
for n= 1: plotPoints + 1
    h= plot(time, realF(:, n), 'LineWidth', 2);
    ylim([-1 1]);
    set(h, 'EraseMode', 'xor');
    frame= getframe(gca);
    aviobj= addframe(aviobj,frame);
end
aviobj= close(aviobj);

Intensity*g*Pressure
Intervals
A_Orders

toc