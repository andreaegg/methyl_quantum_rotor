function [t,signal,ham] = ex_3pESEEM_mqr_1uncoupled_basis_CH3_parallel(Sys,Exp,Opt)

% 3pESEEM simulation script
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
%                    default: 0.2
%
% Exp - struct with fields specifing the experimental parameters
%       .mwfrq       microwave frequency [GHz]
%                    default: 35 (~ Q-Band)
%       .B0          magnetic field [T]
%                    default: 1.224 (= Q-Band)
%       .tau         first interpulse delay [us]
%                    default: 0.120
%       .dt          time increment [us]
%                    default: 0.012
%       .T           starting value for second interpulse delay
%                    default: 0.012
%       .npoints     number of points of time-domain signal
%                    default: 1024
%       .tpi2        pulse length pi/2 pulse [us]
%                    default: 0.012
%       .phase       phase cycle ("+-xy", column = phases, rows = pulses)
%                    default: ["+x","-x"]
%       .addphase    addition of phase cycle signals ("+","-", column = sign, rows = pulses)
%                    default: ["+","-"]
%
% Opt - struct with simulation options
%       .knots       number of orientations on a meridian for with spectrum is simulated
%                    default: 31
%       .weights     weights for orientation selection simulations
%       .vecs        phi and theta dependent vectors for orientation
%                    selection simulations

% Output:
% t          time axis of 3p-ESEEM signal [us]
% signal     ESEEM signal (complex, includes unmodulated part)

% adapted from G. Jeschke and D. Klose 2018
% modified by Andrea Eggeling 2021


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
        Sys.vt = 0.2;
    else
        validateattributes(Sys.vt,{'numeric'},{'nonnegative','scalar'})
    end
end

% if length(find(Sys.Itype == "2H")) >= 3 && isequal(find(Sys.Itype == "2H"),[1 2 3])
%     error('This program simulates the methyl quantum rotor tunnel effect for a CH3-group!')
% end

% Exp input validation %
% -------------------- %

if ~exist('Exp','var')
    Exp.mwfrq   = 35;
    Exp.B0      = 1.224;
    Exp.tau     = 0.120;
    Exp.dt      = 0.012;
    Exp.T       = 0.012;
    Exp.npoints = 1024;
    Exp.tpi2    = 0.012;
    Exp.phase   = ["+x" "+x" "+x" "+x";...
                   "+x" "-x" "+x" "-x";...
                   "+x" "+x" "-x" "-x"];
    Exp.addphase= ["+" "-" "-" "+"];
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
    if ~isfield(Exp,'T')
        Exp.T = 0.012;
    else
        validateattributes(Exp.T,{'numeric'},{'nonnegative','scalar'})
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
    if ~isfield(Exp,'phase')
        Exp.phase = ["+x" "+x" "+x" "+x";...
                     "+x" "-x" "+x" "-x";...
                     "+x" "+x" "-x" "-x"];
    else
        validateattributes(Exp.phase,{'string'},{'nonempty'})
    end
    if ~isfield(Exp,'addphase')
        Exp.addphase = ["+" "-" "-" "+"];
    else
        validateattributes(Exp.addphase,{'string'},{'nonempty'})
    end
    if size(Exp.phase,1) > (3*length(Exp.tpi2))
        error('Number of phase cycled pulses cannot exeed number of applied pulses in pulse sequence.')
    end
    if size(Exp.addphase,1) > (3*length(Exp.tpi2))
        error('Number of phase cycled pulses cannot exeed number of applied pulses in pulse sequence.')
    end
end

% Opt input validation %
% -------------------- %

if ~exist('Opt','var')
    Opt.knots = 31;
else
    if ~isfield(Opt,'knots')
        Opt.knots = 31;
    else
        validateattributes(Opt.knots,{'numeric'},{'nonnegative','scalar'})
    end
    if ~isfield(Opt,'weights')
        Opt.weights = [];
    else
        validateattributes(Opt.weights,{'numeric'},{'nonnegative','vector'})
    end
    if ~isfield(Opt,'vecs')
        Opt.vecs = [];
    else
        validateattributes(Opt.vecs,{'numeric'},{'nrows',3})
    end
end

% Initialize output vectors %
% ------------------------- %

npoints = Exp.npoints;
t = linspace(Exp.T,Exp.T + (Exp.npoints-1)*Exp.dt,Exp.npoints);
signal = zeros(1,Exp.npoints);

% Calculate constants, frequencies, distances and anisotropic HF coupling %
% ----------------------------------------------------------------------- %

ye = Sys.g*bmagn/hbar;

if Sys.Inum > 0
    for k = 1:Sys.Inum
        I(k)      = nucdata(Sys.Itype{k});
        gyro      = nucgval(Sys.Itype{k})*nmagn/hbar;
        yi(k)     = gyro/2/pi/1e6;                                         % [MHz T^-1]
        wI(k)     = -gyro*Exp.B0/2/pi/1e6;                                 % [MHz]
        distance  = norm(Sys.Icoord(k,:) - Sys.Scoord);
        r(k)      = distance;                                              % [Å]
        uv_dd(k,:) = (Sys.Icoord(k,:) - Sys.Scoord)/distance;               % unit vector
        wdd(k)    = ye*gyro*mu0*hbar/(4*pi*(distance*1e-10)^3)/(2*pi*1e6); % [MHz]
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
        rx   = sop(1,'x');
        ry   = sop(1,'y');
        e1   = [1 0 0;0 0 0;0 0 0];
        e2   = [0 0 0;0 1 0;0 0 0];
        e3   = [0 0 0;0 0 0;0 0 1];
        ebig = eye(2^4);
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
            ham = ham + wI(k)*iz{k} + a*sziz{k};
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
            ham0 = kron(e1,hamr1) + kron(e2,hamr2) + kron(e3,hamr3);
            ham = ham0 + to;
        end
    end
    
    % Apply pulse sequence including phase cycle %
    % ------------------------------------------ %
    
    csignal = zeros(1,Exp.npoints);
    
    w_pi2 = 1/(4*Exp.tpi2);
    p1_pc = phasecycle(w_pi2,Exp.phase(1,:),{sx,sy});       % generate pi/2 pulse with phase cycle as cell
    p2_pc = phasecycle(w_pi2,Exp.phase(2,:),{sx,sy});       % generate pi/2 pulse with phase cycle as cell
    p3_pc = phasecycle(w_pi2,Exp.phase(3,:),{sx,sy});       % generate pi/2 pulse with phase cycle as cell
   
    % Generation of propagators
    [upi2_1,upi2_1r] = gen_propagators(Exp.tpi2,p1_pc);
    [upi2_2,upi2_2r] = gen_propagators(Exp.tpi2,p2_pc);
    [upi2_3,upi2_3r] = gen_propagators(Exp.tpi2,p3_pc);
    utau  = expm(-1i*2*pi*ham*Exp.tau);
    utaur = expm( 1i*2*pi*ham*Exp.tau);
    
    % Propagation through pulse sequence pi/2 - tau - pi/2 - T(+dt) - pi/2 - tau - det
    
    for p = 1:npoints
        sigma = zeros(size(sig0));
        for q = 1:length(Exp.addphase)
            uT      = expm(-1i*2*pi*ham*(Exp.T+(p-1)*Exp.dt));
            uTr     = expm(1i*2*pi*ham*(Exp.T+(p-1)*Exp.dt));
            sig_tmp = upi2_1{q}*sig0*upi2_1r{q};          % pi/2 pulse
            sig_tmp = utau*sig_tmp*utaur;                 % tau1
            sig_tmp = upi2_2{q}*sig_tmp*upi2_2r{q};       % pi/2 pulse        
            sig_tmp = uT*sig_tmp*uTr;                     % T initial
            sig_tmp = upi2_3{q}*sig_tmp*upi2_3r{q};       % pi/2 pulse  
            sig_tmp = utau*sig_tmp*utaur;                 % tau2
            switch Exp.addphase(q)
                case "+"
                    sigma = sigma + sig_tmp/length(Exp.addphase);
                case "-"
                    sigma = sigma - sig_tmp/length(Exp.addphase);
            end
        end
        
        csignal(p) = trace(sm*sigma);                   % det
    end
    signal  = signal + weights(ori)*csignal;
end
signal   = signal/sum(weights);
end