function [t,signal,frq,spectrum] = ex_2pESEEM_mqr_parallel(Sys,Exp,Opt)



% 2pESEEM simulation script
%
% takes the following interactions into account:
% - Electron Zeemann
% - Nuclear Zeemann
% - Electron-Nuclear Hyperfine coupling (isotropic and anisotropic)
% - Methyl quantum rotor mixing
% -> negligible dipolar proton-proton coupling
% -> ZFS not taken into account

% requires EasySpin

% Input:
% Sys - struct with fields specifing spin system
%       .g           g-value of observed e-spin
%                    default: 2.000
%       .Scoord      coordinates of observed e-spin [x y z]
%                    default: [0 0 0]
%       .ws          e-spin resonance off-set [MHz]
%                    default: 0
%       .Inum        number of nuclear spin(s)
%                    default = 0
%       .Itype       I nucleus e.g. "1H","14N" in a cell
%                    --> if MQR should be considered protons must be written in positions 1,2,3
%       .Icoord      coordinates of nuclear spin(s) (rows = nucleus, columns = [x y z])
%                    --> if MQR should be considered protons must be written in positions 1,2,3
%       .HFiso       isotropic HF contribution [MHz] (rows = Aiso , column = nucleus)
%                    default = 0
%       .methyl      1 if a methyl group is in close proximity of e-spin
%                    default = 0
%       .vt          methyl quantum rotor tunnel frequency [MHz]
%                    default: 0.1
%
% Exp - struct with fields specifing the experimental parameters
%       .mwfrq       microwave frequency [GHz]
%                    default: 35 (~ Q-Band)
%       .B0          magnetic field [T]
%                    default: 1.224 (= Q-Band)
%       .tau         initial interpulse delay [us]
%                    default: 0.120
%       .dt          time increment [us]
%                    default: 0.012
%       .npoints     number of points of time-domain signal
%                    default: 1024
%       .tpi2        pulse length pi/2 pulse [us]
%                    default: 0.012
%       .tpi         pulse length pi pulse [us]
%                    default: 0.024
%       .phase       phase cycle ("+-xy", column = phases, rows = pulses)
%                    default: ["+x","-x"]
%       .addphase    addition of phase cycle signals ("+","-", column = sign, rows = pulses)
%                    default: ["+","-"]
%
% Opt - struct with simulation options
%       .knots       number of orientations on a meridian for with spectrum is simulated
%                    default: 31
%       .zerofilling dimension of zerofilling for DFT
%                    default: Exp.npoints
%       .apodization apodization window for fourier transformation (see dft_ctav.m)
%                    default: 'dolph-chebyshev'

% Output:
% t          time axis of 2p-ESEEM signal [us]
% signal     ESEEM signal (complex, includes unmodulated part)
% frq        frequency axis of spectrum [MHz]
% spectrum   ESEEM spectrum (complex, negative and positive frequencies)

% adapted from G. Jeschke and D. Klose 2018



% Sys input validation %
% -------------------- %

if ~exist('Sys','var')
    Sys.g      = 2;
    Sys.Scoord = [0 0 0];
    Sys.ws     = 0;
    Sys.Inum   = 0;
    Sys.methyl = 0;
else
    if ~isfield(Sys,'g')
        Sys.g = 2;
    else
        validateattributes(Sys.g,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Sys,'Scoord')
        Sys.Scoord = [0 0 0];
    else
        validateattributes(Sys.Scoord,{'numeric'},{'numel',3})
    end
    if ~isfield(Sys,'ws')
        Sys.ws = 0;
    else
        validateattributes(Sys.ws,{'numeric'},{'scalar'})
    end
    if ~isfield(Sys,'Inum')
        Sys.Inum = 0;
    else
        validateattributes(Sys.Inum,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Sys,'methyl')
        Sys.methyl = 0;
    else
        validateattributes(Sys.methyl,{'numeric'},{'nonnegative','scalar','<=',1})
    end
end

if Sys.Inum > 0
    if isfield(Sys,'Itype')
        if length(Sys.Itype) ~= Sys.Inum
            error('# I spins must all get information on I type. \n Length(Sys.Itype) must be equal to Sys.Inum')
        end
    else
        error('Information on I type must be passed by variable Sys.Itype in as a string e.g. 1H,2H,14N')
    end
    if isfield(Sys,'Icoord')
        if size(Sys.Icoord,1) ~= Sys.Inum
            error('I spins must all get coordinates. \n Rows(Sys.Icoord) must be equal to Sys.Inum.')
        elseif size(Sys.Icoord,2) ~= 3
            error('I spins must all get x,y,z coordinates. \n Columns(Sys.Icoord) must be equal to x,y,z.')
        end
    else
        error('Information on I coordinates must be passed by variable Sys.Icoord in as a matrix (rows = nucleus, columns = [x y z])')
    end
    if isfield(Sys,'HFiso')
        if length(Sys.HFiso) ~= Sys.Inum
            error('# I spins must all get information on isotropic HF interaction.')
        end
    end
    if ~isfield(Sys,'HFiso')
        Sys.HFiso = zeros(1,Sys.Inum);
    end
else
    Sys.Itype  = [];
    Sys.Icoord = [];
    Sys.HFiso  = [];
end

if Sys.methyl == 1
    if ~isfield(Sys,'vt')
        Sys.vt = 1.2;
    else
        validateattributes(Sys.vt,{'numeric'},{'nonnegative','scalar'})
    end
end

% Exp input validation %
% -------------------- %

if ~exist('Exp','var')
    Exp.mwfrq   = 35;
    Exp.B0      = 1.224;
    Exp.tau     = 0.120;
    Exp.dt      = 0.012;
    Exp.npoints = 1024;
    Exp.tpi2    = 0.012;
    Exp.tpi     = 0.024;
    Exp.phase   = ["+x" "-x"];
    Exp.addphase= ["+" "-"];
else
    if ~isfield(Exp,'mwfrq')
        Exp.mwfrq = 35;
    else
        validateattributes(Exp.mwfrq,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Exp,'B0')
        if isfield(Exp,'mwfrq')
            Exp.B0 = planck*1e9*Exp.mwfrq/(Sys.g*bmagn);
        else
            Exp.B0 = 1.224;
        end
    else
        validateattributes(Exp.B0,{'numeric'},{'nonnegative','scalar'})
        Exp.mwfrq = Exp.B0*Sys.g*bmagn/(planck*1e9);
    end
    if ~isfield(Exp,'tau')
        Exp.tau = 0.120;
    else
        validateattributes(Exp.tau,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Exp,'dt')
        Exp.dt = 0.012;
    else
        validateattributes(Exp.dt,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Exp,'npoints')
        Exp.npoints = 1024;
    else
        validateattributes(Exp.npoints,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Exp,'tpi2')
        Exp.tpi2 = 0.012;
    else
        validateattributes(Exp.tpi2,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Exp,'tpi')
        Exp.tpi = 0.012;
    else
        validateattributes(Exp.tpi,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Exp,'phase')
        Exp.phase = ["+x" "-x"];
    else
        validateattributes(Exp.phase,{'string'},{'nonempty'})
    end
    if ~isfield(Exp,'addphase')
        Exp.addphase = ["+" "-"];
    else
        validateattributes(Exp.addphase,{'string'},{'nonempty'})
    end
    if size(Exp.phase,1) > (length(Exp.tpi2)+length(Exp.tpi))
        error('Number of phase cycled pulses cannot exeed number of applied pulses in pulse sequence.')
    end
    if size(Exp.addphase,1) > (length(Exp.tpi2)+length(Exp.tpi))
        error('Number of phase cycled pulses cannot exeed number of applied pulses in pulse sequence.')
    end
    
end

% Opt input validation %
% -------------------- %

if ~exist('Opt','var')
    Opt.knots = 31;
    Opt.apodization = 'dolph-chebyshev';
    Opt.zerofilling = Exp.npoints;
else
    if ~isfield(Opt,'knots')
        Opt.knots = 31;
    else
        validateattributes(Opt.knots,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Opt,'apodization')
        Opt.apodization = 'dolph-chebyshev';
    else
        validateattributes(Opt.apodization,{'char'},{'nonempty'})
        allowedInputs = {'none','hamming','dolph-chebyshev','kaiser','lorentz-gauss','sinebell'};
        validatestring(Opt.apodization,allowedInputs);
    end
    if ~isfield(Opt,'zerofilling')
        Opt.zerofilling = Exp.npoints;
    else
        validateattributes(Opt.zerofilling,{'numeric'},{'nonnegative','scalar'})
    end
end


% Initialize output vectors %
% ------------------------- %

npoints = Exp.npoints;
t = linspace(Exp.tau,(Exp.npoints-1)*Exp.dt,Exp.npoints);
signal = zeros(1,Exp.npoints);
if Opt.zerofilling == Exp.npoints
    spectrum = zeros(1,2*Exp.npoints);
else
    spectrum = zeros(1,Opt.zerofilling);
end

% Calculate constants, frequencies, distances and anisotropic HF coupling %
% ----------------------------------------------------------------------- %

ye = gfree*bmagn/hbar;

if Sys.Inum > 0
    for k = 1:Sys.Inum
        I(k)      = nucdata(Sys.Itype{k});
        gyro      = nucgval(Sys.Itype{k})*nmagn/hbar;
        yi(k)     = gyro/2/pi/1e6;                                         % [MHz T^-1]
        wI(k)     = -gyro*Exp.B0/2/pi/1e6;                                 % [MHz]
        distance  = norm(Sys.Icoord(k,:) - Sys.Scoord);
        r(k)      = distance;                                              % [Å]
        uvec(k,:) = (Sys.Icoord(k,:) - Sys.Scoord)/distance;               % unit vector
        wdd(k)    = ye*gyro*mu0*hbar/(4*pi*(distance*1e-10)^3)/(2*pi*1e6); % [MHz]
    end
    fprintf(1,'The S-H distances are : %4.1f, %4.1f, %4.1f A\n',r(1:3))
end

% Generate spin operators %
% ----------------------- %

multiplenuc = (Sys.Inum > 0);
switch multiplenuc
    case 0
        sx = sop(1/2,'x');
        sy = sop(1/2,'y');
        sz = sop(1/2,'z');
        sm = sx - 1i*sy;
        sig0 = -sz;
    case 1
        spinvec   = [1/2 I];
        stringvec = repmat('e',1,Sys.Inum);
        sx   = sop(spinvec,strcat('x',stringvec));
        sy   = sop(spinvec,strcat('y',stringvec));
        sz   = sop(spinvec,strcat('z',stringvec));
        for k = 1:Sys.Inum
            stringvec      = repmat('e',1,Sys.Inum+1); % +1 because electron spin must also be included
            stringvec(k+1) = 'x';
            ix{k} = sop(spinvec,stringvec);
            
            stringvec      = repmat('e',1,Sys.Inum+1); % +1 because electron spin must also be included
            stringvec(k+1) = 'y';
            iy{k} = sop(spinvec,stringvec);
            
            stringvec      = repmat('e',1,Sys.Inum+1); % +1 because electron spin must also be included
            stringvec(k+1) = 'z';
            iz{k} = sop(spinvec,stringvec);
            
            stringvec      = repmat('e',1,Sys.Inum);
            stringvec(k) = 'z';
            sziz{k} = sop(spinvec,strcat('z',stringvec));
            
            stringvec    = repmat('e',1,Sys.Inum);
            stringvec(k) = 'x';
            szix{k} = sop(spinvec,strcat('z',stringvec));
            
        end
        n    = size(sx,1);
        sm   = sx - 1i*sy;
        sig0 = -sz;
end

% Define Electron spin Hamiltonian in Eigenbasis of coupled nuclei before
% methyl rotor admixing
hamstart = Sys.ws*sz;

if Sys.methyl == 1
    if length(find(Sys.Itype == "1H")) >= 3 && isequal(find(Sys.Itype == "1H"),[1 2 3])
        et   = sop(1,'e');
        sx   = kron(et,sx);
        sy   = kron(et,sy);
        sz   = kron(et,sz);
        sm   = kron(et,sm);
        sig0 = -sz;
        to   = zeros(size(sm));
        for k = 1:2*n
            to(k,k+n) = 1/2;
            to(k+n,k) = 1/2;
        end
        for k = 1:n
            to(k,k+2*n) = 1/2;
            to(k+2*n,k) = 1/2;
        end
    else
        error('Methyl quantum rotor effect can only be taken into account when a methyl group is present in structure, i.e. at least 3 1H required.)');
    end
else
    to = [];
end

% Simulation loop over a set of magnetic field orientations %
% --------------------------------------------------------- %
[vecs,weights] = sphgrid('Ci',Opt.knots,'c');
nori = length(weights); % number of orientations

parfor ori = 1:nori
    
    % Prepare Hamiltonian for the different cases %
    % ------------------------------------------- %
    ham = hamstart;
    sub = zeros(size(sziz{1}));
    addhamr2 = zeros(size(sziz{1}));
    addhamr3 = zeros(size(sziz{1}));
    if Sys.Inum > 0
        % calculate couplings
        cvec = vecs(:,ori)';
        for k = 1:Sys.Inum
            unitvec = uvec(k,:);
            ct  = sum(cvec.*unitvec);
            a   = (3*ct^2-1)*wdd(k) + Sys.HFiso(k);
            b   = 3*ct*sqrt(1-ct^2)*wdd(k);
            ham = ham + wI(k)*iz{k} + a*sziz{k} + b*szix{k};
            if Sys.methyl == 1
                helpvec = repmat([1 2 3],1,2);
                if k <= 3
                    sub  = sub + a*sziz{k};
                    addhamr2 = addhamr2 + a*sziz{helpvec(k+2)};
                    addhamr3 = addhamr3 + a*sziz{helpvec(k+1)};
                end
            end
        end
        if Sys.methyl == 1
            hamr1 = ham;                      % methyl rotation Hamiltonian
            hamr2 = ham - sub + addhamr2;     % methyl rotation Hamiltonian
            hamr3 = ham - sub + addhamr3;     % methyl rotation Hamiltonian
            ham0 = zeros(3*n);
            ham0(1:n,1:n) = hamr1;
            ham0(n+1:2*n,n+1:2*n) = hamr2;
            ham0(2*n+1:3*n,2*n+1:3*n) = hamr3;
            ham = ham0 + 2/3*Sys.vt*to;
        end
    end
    
    % Apply pulse sequence including phase cycle %
    % ------------------------------------------ %
    
    csignal = zeros(1,Exp.npoints);
    
    w_pi2 = 1/(4*Exp.tpi2);
    w_pi  = 1/(2*Exp.tpi);
    
    p1_pc = phasecycle(w_pi2,Exp.phase(1,:),{sx,sy});       % generate pi/2 pulse with phase cycle as cell
    p2_pc = phasecycle(w_pi,[],{sx,sy});                    % generate pi pulse with phase cycle as cell
    
    % Generation of propagators
    upi2 = gen_propagators(ham,Exp.tpi2,p1_pc);
    upi  = gen_propagators(ham,Exp.tpi,p2_pc);
    utau = expm(-1i*2*pi*ham*Exp.tau);
    
    % Propagation trough pulse sequence pi/2 - tau(+dt) - pi - tau(+dt) - det
    sig = propagation_pc(sig0,upi2,Exp.addphase);
    sig = utau*sig*utau';
    for p = 1:npoints
        sigtmp = expm(-1i*2*pi*ham*(p-1)*Exp.dt)*sig*expm(-1i*2*pi*ham*(p-1)*Exp.dt)';
        sigtmp = propagation_pc(sigtmp,upi,[]);
        sigtmp = utau*sigtmp*utau';
        sigtmp = expm(-1i*2*pi*ham*(p-1)*Exp.dt)*sigtmp*expm(-1i*2*pi*ham*(p-1)*Exp.dt)';
        csignal(p) = trace(sm*sigtmp);
    end
    signal  = signal + weights(ori)*csignal;
    csignal = csignal - mean(real(csignal)) - 1i*mean(imag(csignal));
    [cfrq,cspectrum] = dft_ctav(t,csignal,Opt);
    if ori == 1
        frq{ori} = cfrq;
    end
    spectrum = spectrum + weights(ori)*cspectrum;
end
frq = frq{1};
spectrum = spectrum/sum(weights);
signal   = signal/sum(weights);

end