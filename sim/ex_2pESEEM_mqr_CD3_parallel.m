function [t,signal,hamtot] = ex_2pESEEM_mqr_CD3_parallel(Sys,Exp,Opt)

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
%                    default: 0.02
%       .vNQ         nuclear quadrupole coupling frequency for 2H [MHz]
%                    default: 0.2
%       .eta         asymmetry of nuclear quadrupole interaction (0 <= eta <= 1)
%                    default: 0
%       .NQcoord     coordinates used to calculate z(PAS) for nuclear
%                    quadrupolar interaction (rows = nucleus, colums = [x y z] )
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

% Output:
% t          time axis of 2p-ESEEM signal [us]
% signal     ESEEM signal (complex, includes unmodulated part)

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
        Sys.vt = 0.02;
    else
        validateattributes(Sys.vt,{'numeric'},{'nonnegative','scalar'})
    end
end

if numel((Sys.Itype == "2H")) > 0
    if ~isfield(Sys,'vNQ')
        Sys.vNQ = 0.2;
    else
        validateattributes(Sys.vNQ,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Sys,'eta')
        Sys.eta = 0;
    else
        validateattributes(Sys.eta,{'numeric'},{'nonnegative','scalar','<=',1})
    end
    if isfield(Sys,'NQcoord')
        if size(Sys.NQcoord,2) ~= 3
            error('X of X-D bond must all get x,y,z coordinates. \n Columns(Sys.NQcoord) must be equal to x,y,z.')
        end
    else
        error('Deuteriums are in the spin system, therefore the nuclear quadrupole interaction must be taken into account.\n Pass information on z-PAS-axis of X-D.')
    end
end

if length(find(Sys.Itype == "1H")) >= 3 && isequal(find(Sys.Itype == "1H"),[1 2 3])
    error('This program simulates the methyl quantum rotor tunnel effect for a CD3-group!')
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
        Exp.tpi = 0.024;
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

% Initialize output vectors %
% ------------------------- %

npoints = Exp.npoints;
t = linspace(2*Exp.tau,2*Exp.tau +2*(Exp.npoints-1)*Exp.dt,Exp.npoints);
signal = zeros(1,Exp.npoints);

% Calculate constants, frequencies, distances, anisotropic HF coupling and NQI %
% ---------------------------------------------------------------------------- %

ye = gfree*bmagn/hbar;

if Sys.Inum > 0
    for k = 1:Sys.Inum
        I(k)      = nucdata(Sys.Itype{k});
        gyro      = nucgval(Sys.Itype{k})*nmagn/hbar;
        yi(k)     = gyro/2/pi/1e6;                                         % [MHz T^-1]
        wI(k)     = -gyro*Exp.B0/2/pi/1e6;                                 % [MHz]
        distance  = norm(Sys.Icoord(k,:) - Sys.Scoord);
        r(k)      = distance;                                              % [Å]
        uv_dd(k,:)= (Sys.Icoord(k,:) - Sys.Scoord)/distance;               % unit vector
        wdd(k)    = ye*gyro*mu0*hbar/(4*pi*(distance*1e-10)^3)/(2*pi*1e6); % [MHz]
        if numel((Sys.Itype == "2H")) > 0
            K(k) = Sys.vNQ/(4*nucdata(Sys.Itype{k})*(2*nucdata(Sys.Itype{k})-1));   % NQI coupling constant [MHz]
            dist = norm(Sys.Icoord(k,:) - Sys.NQcoord);                             % [Å]
            z    = (Sys.Icoord(k,:) - Sys.NQcoord)/dist;
            xy   = null(z)';
            x    = xy(1,:);
            uv_nqz(k,:) = z;                                                        % z-PAS unit vector
            uv_nqx(k,:) = x;                                                        % x-PAS unit vector
        end
    end
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
            stringvec(k)   = 'z';
            sziz{k} = sop(spinvec,strcat('z',stringvec));
            
            stringvec      = repmat('e',1,Sys.Inum);
            stringvec(k)   = 'x';
            szix{k} = sop(spinvec,strcat('z',stringvec));
            
            inq{k} = 3*iz{k}*iz{k} - (ix{k}*ix{k} + iy{k}*iy{k} + iz{k}*iz{k});
            
        end
        n    = size(sx,1);
        sm   = sx - 1i*sy;
        sig0 = -sz;
end

% Define Electron spin Hamiltonian in Eigenbasis of coupled nuclei before
% methyl rotor admixing
hamstart = Sys.ws*sz;

if Sys.methyl == 1
    if (length(find(Sys.Itype == "2H")) >= 3 && isequal(find(Sys.Itype == "2H"),[1 2 3]))
        et   = sop(1,'e');
        rx   = sop(1,'x');
        ry   = sop(1,'y');
        e1   = [1 0 0;0 0 0;0 0 0];
        e2   = [0 0 0;0 1 0;0 0 0];
        e3   = [0 0 0;0 0 0;0 0 1];
        ebig = eye(2*3^3);
        sx   = kron(et,sx);
        sy   = kron(et,sy);
        sz   = kron(et,sz);
        sm   = kron(et,sm);
        sig0 = -sz;
        rmqr = -(Sys.vt/3)*(rx*rx - ry*ry + sqrt(2)*rx);
        to   = kron(rmqr,ebig);
    else
        error('Methyl quantum rotor effect can only be taken into account when a methyl group is present in structure, i.e. at least 3 1H required.)');
    end
else
    to = [];
end

% Simulation loop over a set of magnetic field orientations %
% --------------------------------------------------------- %
if isempty(Opt.weights)
    [vecs,weights] = sphgrid('Ci',Opt.knots,'c');
elseif ~isempty(Opt.weights)
    weights = Opt.weights;
    vecs    = Opt.vecs;
end

nori = length(weights); % number of orientations

for ori = 1:nori
    P = 0;
    % Prepare Hamiltonian for the different cases %
    % ------------------------------------------- %
    ham = hamstart;
    if Sys.Inum > 0
        sub = zeros(size(sziz{1}));
        addhamr2 = zeros(size(sziz{1}));
        addhamr3 = zeros(size(sziz{1}));
        % calculate couplings
        cvec = vecs(:,ori)';
        for k = 1:Sys.Inum
            uv_curr = uv_dd(k,:);
            ct_dd  = sum(cvec.*uv_curr);
            a   = (3*ct_dd^2-1)*wdd(k) + Sys.HFiso(k);
            % b   = 3*ct_dd*sqrt(1-ct_dd^2)*wdd(k);
            ham = ham + wI(k)*iz{k} + a*sziz{k};% + b*szix{k};
            if Sys.Itype(k) == "2H"
                uv_z = uv_nqz(k,:);
                uv_x = uv_nqx(k,:);
                ct_theta = sum(cvec.*uv_z);
                proj = cvec - uv_x;
                ct_phi = sum(cvec.*proj)/norm(proj);
                P = (3*ct_theta^2-1) + Sys.eta*(1-ct_theta^2)*(ct_phi^2 - (1-ct_phi^2));
                P = K(k)*P;
                ham = ham + P*inq{k};
            end
            if Sys.methyl == 1
                helpvec = repmat([1 2 3],1,2);
                if k <= 3
                    sub  = sub + a*sziz{k} + P*inq{k};% + b*szix{k}
                    addhamr2 = addhamr2 + a*sziz{helpvec(k+2)} + P*inq{helpvec(k+2)};% + b*szix{helpvec(k+2)}
                    addhamr3 = addhamr3 + a*sziz{helpvec(k+1)} + P*inq{helpvec(k+1)};% + b*szix{helpvec(k+1)} 
                end
            end
        end
        if Sys.methyl == 1
            hamr1 = ham;                      % methyl rotation Hamiltonian 0°
            hamr2 = ham - sub + addhamr2;     % methyl rotation Hamiltonian 120°
            hamr3 = ham - sub + addhamr3;     % methyl rotation Hamiltonian 240°
            ham0 = kron(e1,hamr1) + kron(e2,hamr2) + kron(e3,hamr3);
            ham = ham0 + to;
            hamtot = hamr1 + hamr2 + hamr3;
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
    [upi2,upi2r] = gen_propagators(Exp.tpi2,p1_pc);
    [upi,upir]  = gen_propagators(Exp.tpi,p2_pc);
    utau = expm(-1i*2*pi*ham*Exp.tau);
    utaur = expm(1i*2*pi*ham*Exp.tau);
    udt  = expm(-1i*2*pi*ham*Exp.dt);
    udtr = expm(1i*2*pi*ham*Exp.dt);
    
    % Propagation trough pulse sequence pi/2 - tau(+dt) - pi - tau(+dt) - det
    sig = propagation_pc(sig0,upi2,upi2r,Exp.addphase);                 % pi/2 pulse
    sig = utau*sig*utaur;                                               % tau1
    for p = 1:npoints
        sig_tmp = (udt^(p-1))*sig*(udtr^(p-1))  ;                       % +dt (tau1)
        sig_tmp = propagation_pc(sig_tmp,upi,upir,[]);                  % pi pulse
        sig_tmp = utau*sig_tmp*utaur;                                   % tau2
        sig_tmp = (udt^(p-1))*sig_tmp*(udtr^(p-1));                     % + dt (tau2)
        csignal(p) = trace(sm*sig_tmp);                                 % det
    end
    signal  = signal + weights(ori)*csignal;
end
signal   = signal/sum(weights);
end