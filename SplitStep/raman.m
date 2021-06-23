% written by Einat Frishman for Quantum Molecular Dynamics, Weizmann, 2000.
% This program accompanies Ch. 14 of Introduction to Quantum Mechanics:
% A Time Dependent Pespective, by David Tannor (University Science Books, 2001).

% ====== RAMAN

clear
 
% A program to calculate the resonance Raman spectrum for a harmonic ground state
% potential and a displace harmonic excited state potential.
% The evolution of the wave packet is calculated using the split operator propagator.
%
% parameters of the calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%   MOLECULAR PARAMETERS

                          % ===== changeable parameters ============
%mu=0.5*1836;             % the reduced mass
mu=5;
omega=1;
                          % ===== don't change below this line =====
alpha=mu*omega;           %

%%%%%%%%%%%%%%%%%%%%%%%%%%%   X DIMENSION

                          % ===== changeable parameters ============
nop=32;                   % no of points /2
x0=3.0/sqrt(alpha);       % shift of initial wavefunction
xdim=12.0/sqrt(alpha);    % range of x dimension
                          % ===== don't change below this line =====
xrange=nop/xdim;          %
dx=xdim/nop;              %
x=[-1*nop:nop-1]/xrange;  % x range vector
n=length(x);              %
x1=x+x0;                  %
sx1=sqrt(alpha)*(x+x0);   %
delta=sqrt(alpha)*x0;     %

%%%%%%%%%%%%%%%%%%%%%%%%%%%   T DIMENSION

                          % ===== changeable parameters ============
t_total=14*pi/omega;      % total propagation time
                          % a long enough period is required to
                          % calculate the spectrum accurately.
                          % spacing between adjacent k values
                          % is inversly proportional to this quantity.
t_display=4*pi/omega;     % maximum time shown in plots of C_0j.
niter=250;                % no of iterations
nstep=40 ;                % no of time steps per iteration

                          % ===== don't change below this line =====
dt1=t_total/(niter-1);    % time step of each iteration
dt=dt1/nstep;             % time step for propagation
c=-i;                     % constant copmlex -i
t=[0:niter-1]*dt1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shifted gaussian - for ground PES centered at x=-x0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pse=exp(-0.5*sx1.*sx1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   first ten Hermite polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        h0=ones(n,1)';
        h1=2.*sx1;
        h2=(4.*sx1.*sx1-2);
        h3=2*sx1.*h2-2*2.*h1;
        h4=2*sx1.*h3-2*3.*h2;
        h5=2*sx1.*h4-2*4.*h3;
        h6=2*sx1.*h5-2*5.*h4;
        h7=2*sx1.*h6-2*6.*h5;
        h8=2*sx1.*h7-2*7.*h6;
        h9=2*sx1.*h8-2*8.*h7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wave functions associated with the first ten levels 
%             of harmonic oscillator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h0=h0.*pse; hsqw=h0.*conj(h0); normh=(hsqw*ones(n,1)*dx); h0=h0/sqrt(normh);
h1=h1.*pse; hsqw=h1.*conj(h1); normh=(hsqw*ones(n,1)*dx); h1=h1/sqrt(normh);
h2=h2.*pse; hsqw=h2.*conj(h2); normh=(hsqw*ones(n,1)*dx); h2=h2/sqrt(normh);
h3=h3.*pse; hsqw=h3.*conj(h3); normh=(hsqw*ones(n,1)*dx); h3=h3/sqrt(normh);
h4=h4.*pse; hsqw=h4.*conj(h4); normh=(hsqw*ones(n,1)*dx); h4=h4/sqrt(normh);
h5=h5.*pse; hsqw=h5.*conj(h5); normh=(hsqw*ones(n,1)*dx); h5=h5/sqrt(normh);
h6=h6.*pse; hsqw=h6.*conj(h6); normh=(hsqw*ones(n,1)*dx); h6=h6/sqrt(normh);
h7=h7.*pse; hsqw=h7.*conj(h7); normh=(hsqw*ones(n,1)*dx); h7=h7/sqrt(normh);
h8=h8.*pse; hsqw=h8.*conj(h8); normh=(hsqw*ones(n,1)*dx); h8=h8/sqrt(normh);
h9=h9.*pse; hsqw=h9.*conj(h9); normh=(hsqw*ones(n,1)*dx); h9=h9/sqrt(normh);

%norm=(alpha/pi)^0.25;
%h0=h0.*pse/sqrt(dx*2^0*factorial(0))*norm;
%h1=h1.*pse/sqrt(dx*2^1*factorial(1))*norm;
%h2=h2.*pse/sqrt(dx*2^2*factorial(2))*norm;
%h3=h3.*pse/sqrt(dx*2^3*factorial(3))*norm;
%h4=h4.*pse/sqrt(dx*2^4*factorial(4))*norm;
%h5=h5.*pse/sqrt(dx*2^5*factorial(5))*norm;
%h6=h6.*pse/sqrt(dx*2^6*factorial(6))*norm;
%h7=h7.*pse/sqrt(dx*2^7*factorial(7))*norm;
%h8=h8.*pse/sqrt(dx*2^8*factorial(8))*norm;
%h9=h9.*pse/sqrt(dx*2^9*factorial(9))*norm;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hij=[ h0; h1; h2; h3; h4; h5; h6; h7; h8; h9];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);clf;
for j=1:10
  subplot(3,4,j)
  plot(sqrt(alpha)*x,hij(j,1:n)) 
  set(gca,'XLim',[sqrt(alpha)*x(1) sqrt(alpha)*x(end)])
  hold on
  plot(-sqrt(alpha)*x0,0,'k*');
  plot(+sqrt(alpha)*x0,0,'k*');
  drawnow
  xlim=get(gca,'XLim');ylim=get(gca,'YLim');
  str=['\psi_{' num2str(j-1) '}(X)'];
  xx=xlim(1)+0.60*(xlim(2)-xlim(1));yy=ylim(1)+0.85*(ylim(2)-ylim(1));
  text(xx,yy,str);
  if (j==9 | j==10);xlabel('\alpha^{1/2} X');end;
end
%text(-10,6.7,'Eigenfunctions of the Harmonic Oscillator','fontsize',20)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate k vector and vectors of kinetic and potential energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        ak1=[0:nop-1]*pi/xdim;
        ak2=[-1*nop:-1]*pi/xdim;
        ak=[ak1 ak2];
        kine=1.*ak.*ak/(2.*mu);
        ekine=exp(c*dt.*kine);
        pote=0.5*mu*omega^2.*x.*x;     % excited states PES centered at x=0
        epote=exp(c*dt.*pote);
        epote2=exp(.5*c*dt.*pote);
        epotei2=exp(-.5*c*dt.*pote);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial function for propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psinit=(alpha/pi)^(1/4)*exp(-0.5*alpha.*(x+x0).*(x+x0));
psisqr=psinit.*conj(psinit); norm=(psisqr*ones(n,1)*dx);
psinit=psinit/sqrt(norm);
psi=psinit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial correlation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:10
   hram(j,1:n)=hij(j,1:n).*conj(psinit);
   hram1(j)=hram(j,1:n)*ones(n,1);
   hraman(j,1)=hram1(j)*dx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate psi(t) for iter=1,...,niter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	h=figure(2);
    set(h,'position',[0, 50,500,460])%----------------------------
    clf;

        for iter=2:niter         % niter "big" steps - i.e. output points
          psi=psi.*epote2;       % first 1/2 step with potential
          for istep=1:nstep      % nstep "small" steps (alternating kine/pote)
  	    psi=ifft(psi);
	    psi=psi.*ekine;
	    psi=fft(psi);
	    psi=psi.*epote;
          end

	  psi=psi.*epotei2;      % 1/2 step backwards with potential
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate C_0j(t) for iter=1,...,niter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  for j=1:10
	    hram(j,1:n)=hij(j,1:n).*conj(psi(1:n));
	    hram1(j)=hram(j,1:n)*ones(n,1);
	    hraman(j,iter)=hram1(j)*dx;
	  end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot psi(x,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y0=max(h0)*1.1;
hold off;
plot(-sqrt(alpha)*x0,0,'k*');
hold on;
plot(+sqrt(alpha)*x0,0,'k*');
plot(sqrt(alpha)*x,real(psi));
set(gca,'YLim',[-y0,y0]);
xlabel('\alpha^{1/2} X','fontsize',13);
ylabel('\psi(X)','fontsize',13);
plot(sqrt(alpha)*x,imag(psi))
plot(sqrt(alpha)*x,abs(psi),'r')   % plotted for iter=1,...,niter.
drawnow

M(iter-1) = getframe(h);%*************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	end       % --- end of loop over niter.

 %   movie2avi(f,'rmn2.avi','compression','cinepak');%=*****************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot C_0j(t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);clf;
for j=1:10;
subplot(3,4,j);
vt=[1:ceil(length(t)*t_display/t_total)];
plot(t(vt),abs(hraman(j,vt)),'b-')
xlim=[t(1) t(vt(end))];ylim=get(gca,'YLim');
set(gca,'XLim',xlim);
str=['|C_{0' num2str(j-1) '}(t)|'];
xx=xlim(1)+0.05*(xlim(2)-xlim(1));yy=ylim(1)+0.85*(ylim(2)-ylim(1));
text(xx,yy,str);
if (j==9 | j==10);xlabel('t');end;
end;
%text(0,1.5,'Correlation Function - Numeric');
%print -djpeg -r60 'rmn3.jpeg'; %************************************************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npts=niter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discrete k vector
% minimum value is ~2*pi/[ total propagation time ]
% maximum value is ~2*pi/dt1a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bk=[0:niter-1]*2*pi/(niter*dt1);
hpfs=zeros(npts,1)';
bkmax=bk(end)
width=0.5*omega^2*delta^2
nbkmin=1;
nbkmax=min(niter,ceil(niter*2.0*width/bkmax));
                                   % show only relevant energies.
                                   % maximum possible number is niter.

%   FFT(X,N) is the N-point FFT, padded with zeros if X has less
%   than N points and truncated if it has more.
%   For length N input vector x, the DFT is a length N vector X,
%   with elements
%                    N
%      X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
%                   n=1
%
%   t=[0:niter-1]*dt1;                   % defined previously
%   bk=[0:niter-1]*2*pi/(niter*dt1);     % defined previously
%   hft = \sum_t C_0j(t) exp(-i*bk*t)    % the expression we want

 figure(4);clf;
 for j=1:10;
 hpf(j,:)=fft(hraman(j,:),npts);
 hpfs=hpfs+(abs(hpf(j,:))).^2;
 subplot(3,4,j)
 plot(bk(nbkmin:nbkmax),(abs(hpf(j,nbkmin:nbkmax)).^2),'b-')
 set(gca,'XLim',[bk(nbkmin) bk(nbkmax)]);
 set(gca,'YTick',get(gca,'YLim'));
 xlim=get(gca,'XLim');ylim=get(gca,'YLim');
 str=['|F_{0' num2str(j-1) '}(\omega)|'];
 xx=xlim(1)+0.05*(xlim(2)-xlim(1));yy=ylim(1)+0.85*(ylim(2)-ylim(1));
 text(xx,yy,str);
 if (j==9 | j==10);xlabel('\omega');end;
 end;
 
 
 figure(5);clf;
 plot(bk(nbkmin:nbkmax),hpfs(nbkmin:nbkmax),'b-')
 title('Spectrum summed over n');

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analytic expression for absolute value of correlation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);clf;
for j=1:10;
j1=j-1;
% t is a vector.
correlation(:,j)=                                     .../
exp(-delta^2/2.*(1.-exp(-i*omega*t'))-i*omega*t'/2).* .../
(1-exp(-i*omega*t')).^j1.*                            .../
(-1)^j1*delta^j1/sqrt(2^j1*factorial(j1));            .../
subplot(3,4,j)
vt=[1:ceil(length(t)*t_display/t_total)];
plot(t(vt),abs(correlation(vt,j)));
xlim=[t(1) t(vt(end))];ylim=get(gca,'YLim');
set(gca,'XLim',xlim);
str=['|C_{0' num2str(j-1) '}(t)|'];
xx=xlim(1)+0.05*(xlim(2)-xlim(1));yy=ylim(1)+0.85*(ylim(2)-ylim(1));
text(xx,yy,str);
if (j==9 | j==10);xlabel('t');end;
end;
%text(0,1.5,'Correlation Function - Analytic');

