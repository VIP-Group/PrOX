% =========================================================================
% -- PRojection Onto the conveX hull (PrOX)
% -- JED in large SIMO Simulator
% -------------------------------------------------------------------------
% -- (c) 2017 Christoph Studer and Oscar Castañeda
% -- e-mail: studer@cornell.edu and oc66@cornell.edu
% -------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our 
% -- paper: 
% -- Oscar Castañeda, Tom Goldstein, and Christoph Studer,
% -- "VLSI Designs for Joint Channel Estimation and Data Detection in Large
% -- SIMO Wireless Systems,"
% -- IEEE Transactions on Circuits and Systems I: Regular Papers,
% -- vol. 65, no. 3, pp. 1120-1132, Mar. 2018.
% =========================================================================

function PrOX_SIMO_JED_sim(varargin)

  % -- set up default/custom parameters
  
  if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
        
    % set default simulation parameters     
    par.runId = 0;        % simulation ID (used to reproduce results)
    par.MR = 16;          % receive antennas 
    par.MT = 1;           % transmit antennas (set not larger than MR!) 
    par.Time = 17;        % time slots (K+1)
    par.TimeDL = 1;       % time slots used for downlink communication
    par.mod = 'BPSK';     % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 1e4;     % number of Monte-Carlo trials (transmissions)
    par.simName = ...     % simulation name (used for saving results)
      ['NCERR_' num2str(par.MR) 'x' num2str(par.MT) '_' ...
      num2str(par.Time), 'TS_', par.mod, '_', ...
      num2str(par.trials), 'Trials'];    
    par.SNRdB_list = ...  % list of SNR [dB] values to be simulated
      -12:2:0;
    par.detector = ...    % define detector(s) to be simulated
      {'MRC-CSIR','MRC','MRC-RT','ML-JED','TASER','PROX','APROX'};
    %NOTE: MRC and MRC_RT only differ in the downlink transmission; MRC_RT
    %      uses the received data to a new channel estimate used for beam-
    %      forming, while MRC uses the same channel estimate used during
    %      the uplink transmission.
    
    % channel parameters
    par.los = 0;          % use line-of-sight channel model?
    par.lospwm = 0;       % use planar wave model for LoS? If not, uses
                          % spherical wave model
    
    % algorithm parameters
    par.PrOX.tmax = 5;    % number of PrOX and APrOX iterations
    par.TASER.tmax = 10;  % number of TASER iterations
    par.TASER.aSc = 0.99; % alpha scale for TASER's step size
                          % use a value of 0.99 for square systems (MR==MT)    
    switch (par.mod)
      case 'BPSK' %FOR BPSK, MR=16 and Time=17,   
        par.PROX.aSc = 1;
        par.PROX.rho = 2^3;
        par.APROX.aSc = 1/32;
        par.APROX.rho = 2^1;
      case 'QPSK' %FOR QPSK, MR=16 and Time=17,
        par.PROX.aSc = 1;
        par.PROX.rho = 2^4;
        par.APROX.aSc = 1/32;
        par.APROX.rho = 2^1;
    end                
    
  else
      
    disp('use custom simulation settings and parameters...')    
    par = varargin{1}; % only argument is par structure
    
  end

  % -- initialization
  
  % use runId random seed (enables reproducibility)
  rng(par.runId); 

  % set up Gray-mapped constellation alphabet (according to IEEE 802.11)
  switch (par.mod)
    case 'BPSK',
      par.symbols = [ -1 1 ];
    case 'QPSK', 
      par.symbols = [ -1-1i,-1+1i, ...
                      +1-1i,+1+1i ];
    case '16QAM',
      par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
                      -1-3i,-1-1i,-1+3i,-1+1i, ...
                      +3-3i,+3-1i,+3+3i,+3+1i, ...
                      +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
      par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                      +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
                         
    case '8PSK',
      par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
                      exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
                      exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
                      exp(1i*2*pi/8*4), exp(1i*2*pi/8*5)];
  end

  % extract average symbol energy
  par.Es = mean(abs(par.symbols).^2); 
  % We will fix the downlink vector to have the same energy par.Es
  par.P = par.Es;
  
  % precompute bit labels
  par.Q = log2(length(par.symbols)); % number of bits per symbol
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

  % track simulation time
  time_elapsed = 0;
  
  % -- start simulation 
  
  % initialize result arrays (detector x SNR)
  res.PER = zeros(length(par.detector),length(par.SNRdB_list)); % vector error rate
  res.SER = zeros(length(par.detector),length(par.SNRdB_list)); % symbol error rate
  res.BER = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate
  %For the downlink
  res.PERDL = zeros(length(par.detector),length(par.SNRdB_list)); % vector error rate
  res.SERDL = zeros(length(par.detector),length(par.SNRdB_list)); % symbol error rate
  res.BERDL = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate
  %MSE for the channel estimate
  res.HMSE = zeros(length(par.detector),length(par.SNRdB_list));


  % trials loop
  tic
  for t=1:par.trials
    par.stop = 0;
    % generate random bit stream (antenna x bit x time slots)
    bits = randi([0 1],par.MT,par.Q,par.Time);
    bits(:,:,1)=0; % nail down first signal
  
    % generate transmit symbol
    for time=1:par.Time
      idx(:,time) = bi2de(bits(:,:,time),'left-msb')+1;
      S(time,:) = par.symbols(idx(:,time));
    end
    
    % generate downlink transmit symbol
    bitsDL = randi([0 1],par.MT,par.Q,par.TimeDL);
    for time=1:par.TimeDL
      idxDL(:,time) = bi2de(bitsDL(:,:,time),'left-msb')+1;
      SDL(:,time) = par.symbols(idxDL(:,time));
    end    
  
    % generate iid Gaussian channel matrix & noise matrices    
    N = sqrt(0.5)*(randn(par.MR,par.Time)+1i*randn(par.MR,par.Time)); 
    NDL = sqrt(0.5)*(randn(par.MT,par.TimeDL)+1i*randn(par.MT,par.TimeDL));
    
    % channel    
    if (par.los)
      [H_swm, H_pwm] = los(par);
      if (par.lospwm)
        H = H_pwm.';
      else
        H = H_swm.';
      end
    else
      H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    end

    % CHEST noise    
    NH = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT)); 
    
    % transmit over noiseless channel (will be used later)
    X = H*S';
  
    % SNR loop
    for k=1:length(par.SNRdB_list)
      
      % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
      N0 = par.MT*par.Es*10^(-par.SNRdB_list(k)/10);      
      
      % transmit data over noisy channel
      Y = X+sqrt(N0)*N;    
      Hest = H + sqrt(N0)*NH; % channel under estimation errors
            
      % algorithm loop      
      for d=1:length(par.detector) 
        switch (par.detector{d}) % select algorithms
          case 'MRC-CSIR', % MRC with perfect CSIR (no CHEST errors)
            [idxhat,bithat,htilde] = MRC(par,H,Y,S,0);
          case 'MRC',      % MRC with imperfect CSIR
            [idxhat,bithat,htilde] = MRC(par,Hest,Y,S,0);
          case 'MRC-RT',   % MRC with imperfect CSIR that retrains channel 
                           % with the estimated uplink data
            [idxhat,bithat,htilde] = MRC(par,Hest,Y,S,1);
          case 'PROX',     % PRojection Onto conveX hull (PrOX)
            [idxhat,bithat,htilde] = PROX(par,Y);
          case 'APROX',    % Approximate PRojection Onto conveX hull (APrOX)
            [idxhat,bithat,htilde] = APROX(par,Y);          
          case 'TASER',    % Triangular Approximate 
                           % SEmidefinite Relaxation (TASER)
            [idxhat,bithat,htilde] = TASER(par,Y);                      
          case 'ML-JED',   % ML detection using sphere decoding
            [idxhat,bithat,htilde] = ML(par,Y);            
          otherwise,
            error('par.detector type not defined.')      
        end

        % -- compute error metrics
        err = (idx~=idxhat);
        res.PER(d,k) = res.PER(d,k) + any(err(:));
        res.SER(d,k) = res.SER(d,k) + sum(err(:))/par.MT/par.Time;    
        res.BER(d,k) = res.BER(d,k) + ... 
                         sum(bits(:)~=bithat(:))/(par.MT*par.Time*par.Q);      
      
        % -- beamforming            
        BSDL = sqrt(par.P)*conj(htilde)*SDL/norm(htilde,2)/sqrt(par.Es);
        XDL = H.'*BSDL;
        YDL = XDL+sqrt(N0)*NDL;
          
        mt=1;
        [~,idxhatDL(mt,:)] = min(abs(YDL.'*ones(1,length(par.symbols))-ones(par.TimeDL,1)*par.symbols).^2,[],2);
        bithatDL(mt,:,:) = par.bits(idxhatDL(mt,:),:)';
          
        res.HMSE(d,k) = res.HMSE(d,k) + norm(htilde - H,2)^2;
          
        errDL = (idxDL~=idxhatDL);
        res.PERDL(d,k) = res.PERDL(d,k) + any(errDL(:));
        res.SERDL(d,k) = res.SERDL(d,k) + sum(errDL(:))/par.MT/par.TimeDL;    
        res.BERDL(d,k) = res.BERDL(d,k) + sum(bitsDL(:)~=bithatDL(:))/(par.MT*par.TimeDL*par.Q);                        
        
      end % algorithm loop
                 
    end % SNR loop    
    
    % keep track of simulation time    
    if toc>10
      time=toc;
      time_elapsed = time_elapsed + time;
      fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
      tic
    end      
  
  end % trials loop

  % normalize results
  res.PER = res.PER/par.trials;
  res.SER = res.SER/par.trials;
  res.BER = res.BER/par.trials;
  res.PERDL = res.PERDL/par.trials;
  res.SERDL = res.SERDL/par.trials;
  res.BERDL = res.BERDL/par.trials;
  res.time_elapsed = time_elapsed;
  
  % -- save final results (par and res structure)
  save([ par.simName '_' num2str(par.runId) ],'par','res');    
    
  % -- show results (generates fairly nice Matlab plot) 
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:','bs-.'};    
  figure('Name','Uplink')
  for d=1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.SER(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.SER(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('symbol error rate (SER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-3 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)

  figure('Name','Downlink')
  for d=1:length(par.detector)
    if d==1
      semilogy(par.SNRdB_list,res.SERDL(d,:),marker_style{d},'LineWidth',2)
      hold on
    else
      semilogy(par.SNRdB_list,res.SERDL(d,:),marker_style{d},'LineWidth',2)
    end
  end
  hold off
  grid on
  xlabel('average SNR per receive antenna [dB]','FontSize',12)
  ylabel('symbol error rate (SER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-3 1])
  legend(par.detector,'FontSize',12)
  set(gca,'FontSize',12)  
 
end

% -- set of detector functions 

%% Maximum Ratio Combining (MRC)
function [idxhat,bithat,htilde] = MRC(par,H,Y,S,RT)
  Z = Y-H*S';
  for time=1:par.Time
    for m=1:par.MT
      hm = H(:,m);
      yhat = Z(:,time)+hm*S(time,m)';
      xhat(m,1) = yhat'*hm/norm(hm,2)^2;
    end 
    [~,idxhat(:,time)] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
    bithat(:,:,time) = par.bits(idxhat(:,time),:);
  end
  if(RT)
    shat = par.symbols(idxhat).';
    htilde = Y*shat/norm(shat,2)^2;
  else
    htilde = H;
  end    
end

%% PRojection Onto conveX hull (PrOX)
function [idxhat,bithat,hest] = PROX(par,Y)

  % -- preprocessing
  G = Y'*Y;
  alpha = par.PROX.aSc*(max(svd(G))+1);
  Gtilde = inv(eye(par.Time)-1/alpha*G);
  % --- find the elements with the greatest absolute value in the diagonal
  %     under the main diagonal to scale Gtilde, so the greatest absolute
  %     value in the scaled matrix is close to one in absolute value
  underDiag = zeros(size(Gtilde,1)-1,1);
  for i=1:length(underDiag)
      underDiag(i) = Gtilde(i+1,i);
  end
  gamma = (2^(ceil(log2(max(abs(real(underDiag)))))+1));
  % --- rescale Gtilde to obtain Ghat
  Ghat = Gtilde/gamma;
  
  rho = par.PROX.rho;
  % -- initialization
  s = G(:,1)/G(1,1)*par.symbols(1);
  
  % -- PrOX main loop
  for ii=1:par.PrOX.tmax
    s = Ghat*s;
    s = rho*s;
    s = min(max(real(s),-1),1) + ...
        (~strcmp(par.mod,'BPSK'))*1i*min(max(imag(s),-1),1);
    s(1) = par.symbols(1);
  end

  % -- compute hard-output estimates of the transmitted signals
  shat = s;  
  mt=1;
  [~,idxhat(mt,:)] = min(abs(shat*ones(1,length(par.symbols))-ones(par.Time,1)*par.symbols).^2,[],2);
  bithat(mt,:,:) = par.bits(idxhat(mt,:),:)';    
  shat = par.symbols(idxhat(mt,:)).';
  
  % -- compute channel estimate
  hest = Y*shat/norm(shat,2)^2;
  
end

%% Approximate PRojection Onto conveX hull (APrOX)
function [idxhat,bithat,hest] = APROX(par,Y)

  % -- preprocessing
  G = Y'*Y;
  alpha = par.APROX.aSc*(max(svd(G))+1);
  Gtilde = eye(par.Time)+1/alpha*G;
  % --- find the elements with the greatest absolute value in the diagonal
  %     under the main diagonal to scale Gtilde, so the greatest absolute
  %     value in the scaled matrix is close to one in absolute value
  underDiag = zeros(size(Gtilde,1)-1,1);
  for i=1:length(underDiag)
      underDiag(i) = Gtilde(i+1,i);
  end
  gamma = (2^(ceil(log2(max(abs(real(underDiag)))))+1));
  % --- rescale Gtilde to obtain Ghat
  Ghat = Gtilde/gamma;
  
  rho = par.APROX.rho;
  % -- initialization
  s = G(:,1)/G(1,1)*par.symbols(1);
  
  % -- PrOX main loop
  for ii=1:par.PrOX.tmax
    s = Ghat*s;
    s = rho*s;
    s = min(max(real(s),-1),1) + ...
        (~strcmp(par.mod,'BPSK'))*1i*min(max(imag(s),-1),1);
    s(1) = par.symbols(1);
  end

  % -- compute hard-output estimates of the transmitted signals
  shat = s;  
  mt=1;
  [~,idxhat(mt,:)] = min(abs(shat*ones(1,length(par.symbols))-ones(par.Time,1)*par.symbols).^2,[],2);
  bithat(mt,:,:) = par.bits(idxhat(mt,:),:)';    
  shat = par.symbols(idxhat(mt,:)).';
  
  % -- compute channel estimate
  hest = Y*shat/norm(shat,2)^2;
  
end

%% detection via Triangular Approximate SEmidefinite Relaxation (TASER)
function [idxhat,bithat,hest] = TASER(par,Y)

  % -- Re-express the problem as in the Massive MU-MIMO case.
  y = Y(:,1)*par.symbols(1);
  H = Y(:,2:end);
      
  switch par.mod
    case 'BPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];
      % -- preprocessing for SDR  
      T = -[HR'*HR HR'*yR ; yR'*HR yR'*yR];           
    case 'QPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];
      % -- preprocessing for SDR 
      T = -[HR'*HR HR'*yR ; yR'*HR yR'*yR];        
    otherwise
      error('modulation type not supported')
  end
  
  % -- preconditioning for SDR
  DInv = diag(1./sqrt(abs(diag(T))));
  Ttilde = DInv*T*DInv;
  stepsize = 1/norm(Ttilde,2);
  
  % -- use standard gradient on non-convex problem  
  gradf = @(L) 2*tril(L*Ttilde);
  proxg = @(L,t) prox_normalizer(L,diag(DInv).^-1);
  
  % Initialize Ltilde  
  Ltilde = diag(diag(DInv).^-1);  
  
  % -- Fast Iterative Soft Thresholding [Beck & Tebouille, 2009]  
  for k = 1:par.TASER.tmax
    Ltilde = proxg(Ltilde-stepsize*gradf(Ltilde)); % compute proxy    
  end  
  
  % -- post processing  
  sest = sign(Ltilde(end,:))';
  shat(1) = par.symbols(1);
  
  switch par.mod
    case 'BPSK'   
      shat(2:par.Time,1) = sest(1:par.Time-1); 
    case 'QPSK'        
      shat(2:par.Time,1) = sest(1:par.Time-1) + 1i*sest(par.Time:end-1);  
    otherwise
      error('modulation type not supported')
  end
  
  % -- compute outputs  
  mt=1;
  [~,idxhat(mt,:)] = min(abs(shat*ones(1,length(par.symbols))-ones(par.Time,1)*par.symbols).^2,[],2);
  bithat(mt,:,:) = par.bits(idxhat(mt,:),:)'; 
  
  shat = par.symbols(idxhat(mt,:)).';
  hest = Y*shat/norm(shat,2)^2;
  
end

% normalize columns of Z to have norm equal to its corresponding scale
function Q = prox_normalizer(Z,scale)
  [N,~] = size(Z);
  Q = Z.*(ones(N,1)*(1./sqrt(sum(abs(Z).^2,1)).*scale'));
end

%% ML detection using sphere decoding
function [idxhat,bithat,hest] = ML(par,Y)

  % -- Re-express the problem as in the Massive MU-MIMO case.
  y = Y(:,1)*conj(par.symbols(1));
  H = Y(:,2:end);  
      
  switch par.mod
    case 'BPSK'
      % -- convert to real domain
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) ; imag(H) ];
      % -- preprocessing
      T = [HR'*HR HR'*yR ; yR'*HR yR'*yR];
      N = par.Time;
    case 'QPSK'
      % -- convert to real domain                  
      yR = [real(y) ; imag(y) ];
      HR = [ real(H) -imag(H) ; imag(H) real(H) ];
      % -- preprocessing
      T = [HR'*HR HR'*yR ; yR'*HR yR'*yR];
      N = par.Time*2-1;        
    otherwise
      error('modulation type not supported')
  end
  
  symbols = [ -1 1 ];
  
  % -- initialization  
  Radius = inf;
  PA = zeros(N,1); % path
  ST = zeros(N,length(symbols)); % stack  

  % -- preprocessing
  R = chol(max(eig(T)+1)*eye(N)-T);
  
  % -- add root node to stack
  Level = N; 
  ST(Level,:) = abs(R(Level,Level)*symbols.').^2;
  
  % -- begin sphere decoder
  while ( Level<=N )          
    % -- find smallest PED in boundary    
    [minPED,idx] = min( ST(Level,:) );
    
    % -- only proceed if list is not empty
    if minPED<inf
      ST(Level,idx) = inf; % mark child as tested        
      NewPath = [ idx ; PA(Level+1:end,1) ]; % new best path
      
      % -- search child
      if ( minPED<Radius )
        % -- valid candidate found
        if ( Level>1 )                  
          % -- expand this best node
          PA(Level:end,1) = NewPath;
          Level = Level-1; % downstep
          DF = R(Level,Level+1:end) * symbols(PA(Level+1:end,1)).';         
          ST(Level,:) = minPED + abs(R(Level,Level)*symbols.'+DF).^2;
        else
          % -- valid leaf found     
          idxML = NewPath;
          bitML = par.bits(idxML.',:);
          % -- update radius (radius reduction)
          Radius = minPED;    
        end
      end      
    else
      % -- no more childs to be checked
      Level=Level+1;      
    end    
  end  

  % -- post processing
  sest = symbols(idxML).';
  shat(1) = par.symbols(1);
  
  switch par.mod
      case 'BPSK'                
        shat(2:par.Time,1) = par.symbols(1)*sest(1:par.Time-1);                  
      case 'QPSK'        
        shat(2:par.Time,1) = real(par.symbols(1))*sest(1:par.Time-1) ...
                             +1i*imag(par.symbols(1))*sest(par.Time:end-1);
  end
  
  % -- compute outputs  
  mt=1;
  [~,idxhat(mt,:)] = min(abs(shat*ones(1,length(par.symbols))-ones(par.Time,1)*par.symbols).^2,[],2);
  bithat(mt,:,:) = par.bits(idxhat(mt,:),:)'; 
    
  hest = Y*shat/norm(shat,2)^2;
  
end

%% Line of sight channel generation.
%  Special thanks to Sven Jacobsson for this function
function [H_swm, H_pwm] = los(par)

    U = par.MT;
    B = par.MR;

    c = 3e8; % speed of light [m/s]
    lambda =  c / 2e9; % carrier wavelength [m]
    delta = .5; % antenna spacing
 
    % distance spread  
    d_spread = 0; % distance spread
    d_max = 150 + d_spread; % max user distance [m]
    d_min = 150 - d_spread; % min user distance [m]
     
    if d_spread > 0
        d_avg = 2/3*(d_max^3-d_min^3)/(d_max^2-d_min^2); % avg user distance [m]
    else
        d_avg = d_max;
    end
    d_UE = sqrt((d_max^2-d_min^2)*rand(U,1) + d_min^2);% UE dist [m]
    aod_UE = unifrnd(45-60,45+60,U,1); % AoD [deg]
 
    broadside_BS_deg = 45; % broadside of BS antenna array  [deg]
 
    coord_BS = [0, 0]; % BS coord.
    coord_UE = ones(U,1)*coord_BS + (d_UE*ones(1,2)).*[cosd(aod_UE), sind(aod_UE)]; % UE coord.
     
    d_ant_BS = delta * lambda; % distance between BS antenna elements [m]
    d_array_BS = (B - 1) * d_ant_BS; % size of BS antenna array [m]
 
    % array rotation
    Omega_BS_deg = wrapTo360(90 - broadside_BS_deg); % BS array rotation [deg]
    Omega_BS_rad = pi/180 * Omega_BS_deg; % BS array rotation [rad]
    Omega_UE_deg = wrapTo360(unifrnd(-180, 180, U, 1) - 90); % UE array rotation [deg]
    Omega_UE_rad = pi/180 * Omega_UE_deg; % UE array rotation [rad]
 
    % antenna elem. coordinates
    x_BS = coord_BS(1) + d_ant_BS*((1:B) - (B+1)/2)*cos(pi-Omega_BS_rad);
    y_BS = coord_BS(2) + d_ant_BS*((1:B) - (B+1)/2)*sin(pi-Omega_BS_rad);
    x_UE = coord_UE(:,1);
    y_UE = coord_UE(:,2);
 
    % LoS AoD
    theta_LoS_BS_rad = Omega_BS_rad - pi/2 + atan2((coord_UE(:,2) - coord_BS(2)),(coord_UE(:,1) - coord_BS(1)));
    theta_LoS_BS_deg = 180/pi * theta_LoS_BS_rad; 
    theta_LoS_UE = pi - Omega_UE_rad - atan2((coord_UE(:,1) - coord_BS(1)),(coord_UE(:,2) - coord_BS(2)));
    theta_LoS_UE_deg = 180/pi * theta_LoS_UE;
 
    % coordinates
    xx_BS = ones(U,1)*x_BS; yy_BS = ones(U,1)*y_BS;
    xx_UE = x_UE*ones(1,B); yy_UE = y_UE*ones(1,B);
     
    % reference distance
    d_ref = sqrt((xx_BS - xx_UE).^2 + (yy_BS - yy_UE).^2);
 
    % angles
    theta_BS = Omega_BS_rad - pi/2 + atan2((yy_UE - yy_BS),(xx_UE - xx_BS));
    theta_UE = pi - Omega_UE_rad*ones(1,B) - atan2((xx_UE - xx_BS),(yy_UE - yy_BS));
     
    % distances between UE and BS antenna elements
    dd_ant_BS = d_ant_BS*ones(U,1)*((1:B)-1); tt_BS = theta_BS(:,1)*ones(1,B);
     
    % distance according to PWM model
    d_pwm = d_ref(:,1)*ones(1,B) - dd_ant_BS.*sin(tt_BS);
     
    % distance according to SWM model
    d_swm = sqrt((d_ref(:,1)*ones(1,B) - dd_ant_BS.*sin(tt_BS)).^2 + (dd_ant_BS.*cos(tt_BS)).^2);
 
    if sum(sum(abs(d_ref - d_swm) < 0.001./d_ref)) ~= B*U
       warning('something is wrong here...'); 
    end
 
    % channel matrix
    H_pwm = d_avg./d_pwm .* exp(-1i*2*pi*d_pwm/lambda);
    H_swm = d_avg./d_swm .* exp(-1i*2*pi*d_swm/lambda);
     
 
end
 
function lon = wrapTo360(lon)
 
    positiveInput = (lon > 0);
    lon = mod(lon, 360);
    lon((lon == 0) & positiveInput) = 360;
 
end
