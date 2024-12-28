clearvars; clc;
close all;
NNorPDA = 1; % 1 for NN, 2 for PDA

T = 1;
t = 100;
time_steps = 1:T:t;

xrange = 10e3;
yrange = 10e3;
vmax = 50;
kappa = 6;

A = [eye(2),T*eye(2);0*eye(2),eye(2)];
G = [T^2/2*eye(2);T*eye(2)];
C = [1,0,0,0;0,1,0,0];

gate_threshold = chi2inv(0.99,2);
track_deletion_treshold = 3;

N1 = 2;
M2 = 2;
N2 = 3;

x0 = [5000, 5000, 25, 25]';

mu_ak = [0,0]';
Q = [2^2,0;0,2^2];

trueTarget = zeros(4,t/T); 

for i = 1:length(time_steps)
    if i==1
        trueTarget(:,i) = x0;
    else
    ak = mvnrnd(mu_ak,Q)';
    trueTarget(:,i) = A*trueTarget(:,i-1) +G*ak;
    end
end

x_position_true = trueTarget(1,:);
y_position_true = trueTarget(2,:);
x_velocity_true = trueTarget(3,:);
y_velocity_true = trueTarget(4,:);

figure;
subplot(2,2,1);
plot(time_steps,x_position_true,LineWidth=1.5,Color="#A2142F");
title("True x-position of the Target vs. Time");
xlabel("k");
xticks(0:10:100);
ylabel("x position");
ylim([0,xrange]);
legend("x position");
grid on;

subplot(2,2,3);
plot(time_steps,y_position_true,LineWidth=1.5,Color="#D95319");
title("True y-position of the Target vs. Time");

xlabel("k");
xticks(0:10:100);
ylabel("y position");
ylim([0,yrange]);
grid on;
legend("y position");

subplot(2,2,[2,4]);
plot(x_position_true,y_position_true,LineWidth=1.5,Color="#77AC30");
title("True Target Trajectory");
ylabel("y position");
ylim([0,yrange]);
xlabel("x-position");
xlim([0,xrange]);
legend("True Target Trajectory");
grid on;

figure;
subplot(2,1,1);
plot(time_steps,x_velocity_true,LineWidth=1.5,Color="#A2142F");
title("True x-velocity of the Target vs. Time");
xlabel("k");
xticks(0:10:100);
ylabel("x velocity");
legend("x velocity");
grid on;

subplot(2,1,2);
plot(time_steps,y_velocity_true,LineWidth=1.5,Color="#D95319");
title("True y-velocity of the Target vs. Time");
xlabel("k");
xticks(0:10:100);
ylabel("y velocity");
legend("y velocity");
grid on;

beta_FA = 1e-7;
V = xrange * yrange;
false_alarms = cell(2,200);

for i= 1:200
    m_k = poissrnd(beta_FA*V);
    FA_xvalues = xrange*rand(1,m_k);
    FA_yvalues = yrange*rand(1,m_k);
    false_alarms{1,i} = FA_xvalues;
    false_alarms{2,i} = FA_yvalues;
end

measurement_sigma = 20;
P_d = 0.9;
R = measurement_sigma^2*eye(2);
measurement_mu = [0;0];

measurements = zeros(2,length(time_steps));

for i = 1:length(time_steps)
    u = rand();
    if u<=P_d
        %detection
        v_k = mvnrnd(measurement_mu,R)';
        measurement = trueTarget(1:2,i) +v_k;
        measurements(:,i) = measurement;
    else
        %miss
        measurements(1,i) = inf;
        measurements(2,i) = inf;
    end
end

track_start_time = randi([0, 50])+1;

detections = cell(2,200);

for i = 1:200

    if i>=track_start_time && i < track_start_time+length(measurements(1,:))
        measurement = measurements(:,i-track_start_time+1);
        if measurement(1) ~= inf
            detections{1,i} = [false_alarms{1,i} measurement(1)];
            detections{2,i} = [false_alarms{2,i} measurement(2)];
        else
  
            detections{1,i} = false_alarms{1,i};
            detections{2,i} = false_alarms{2,i};
        end
    else 
        detections{1,i} = false_alarms{1,i};
        detections{2,i} = false_alarms{2,i};
    end
end

figure;
hold on;


ConfirmedTracks = [];
AllTracks = [];
%for each time step
for i = 1:200
    current_detections_x = detections{1,i};
    current_detections_y = detections{2,i};
    AssociatedMeasurementSet = {};

    if i == 1
        for k = 1:length(current_detections_x)
            fresh_initiator = Initiator([current_detections_x(k);current_detections_y(k)]...
                ,measurement_sigma^2,vmax,kappa);
            InitiatorArray(k) = fresh_initiator;
        end
    else

        MeasurementSet = cell(1,length(current_detections_x));
        for k = 1:length(current_detections_x)
            MeasurementSet{k} = [current_detections_x(k);current_detections_y(k)];
        end
        OriginalMeasurementSet = MeasurementSet;
        indexes_to_delete = [];
        for k =1:length(ConfirmedTracks)
            [ConfirmedTracks(k),MeasurementSet] = ...
                ConfirmedTracks(k).Gating(MeasurementSet,gate_threshold);
            if NNorPDA == 1
            [ConfirmedTracks(k),associated_measurement] = ConfirmedTracks(k).NN();
            if associated_measurement ~= inf
                AssociatedMeasurementSet{end+1} = associated_measurement;
                % measurement exists
                ConfirmedTracks(k) = ConfirmedTracks(k).measurementUpdate(associated_measurement);
                ConfirmedTracks(k).NumberOfConsecutiveMissedDetections = 0;
                ConfirmedTracks(k) = ConfirmedTracks(k).predictionUpdate();
            else
                % no measurement
                ConfirmedTracks(k) = ConfirmedTracks(k).predictionUpdate();
                ConfirmedTracks(k).NumberOfConsecutiveMissedDetections...
                    = ConfirmedTracks(k).NumberOfConsecutiveMissedDetections + 1;
            end
            elseif NNorPDA == 2
                ConfirmedTracks(k) = ConfirmedTracks(k).PDA(P_d,1.0,beta_FA);
                ConfirmedTracks(k) = ConfirmedTracks(k).predictionUpdate();
            end

            if ConfirmedTracks(k).NumberOfConsecutiveMissedDetections == track_deletion_treshold
                indexes_to_delete = [indexes_to_delete,k];
            end
        end
        ConfirmedTracksAfterDeletion = [];
        for k = 1:length(ConfirmedTracks)
            if ~ismember(k,indexes_to_delete)
                ConfirmedTracksAfterDeletion = [ConfirmedTracksAfterDeletion,ConfirmedTracks(k)];
            else
                AllTracks = [AllTracks,ConfirmedTracks(k)];
            end
        end
        ConfirmedTracks = ConfirmedTracksAfterDeletion; 
        
        GateDecisions = Gating(InitiatorArray,MeasurementSet,T,vmax,gate_threshold,C,R);
        AssociationDecisions = Associate(InitiatorArray,MeasurementSet,GateDecisions);
        InitiatorArray = Update(InitiatorArray,MeasurementSet,...
            AssociationDecisions,A,G,Q,T,measurement_sigma^2,C,R,vmax,kappa);
        [InitiatorArray,NewTracks] = CheckInitiators(InitiatorArray,N1,M2,N2);
        for k = 1:length(NewTracks)
            ConfirmedTracks = [ConfirmedTracks, ...
                ConfirmedTrack(NewTracks(k).CovarianceEstimate,...
                NewTracks(k).StateEstimate,A,G,Q,C,R)];
        end
    end

    scatter(detections{1,i}, detections{2,i},'black');
    for l = 1:length(AssociatedMeasurementSet)
        current_point = AssociatedMeasurementSet{l};
        scatter(current_point(1),current_point(2),'blue')
    end

    title(['Time step: ', num2str(i)]);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    xlim([0, xrange]);
    ylim([0, yrange]); 

    % Update the plot
    drawnow;


    pause(0.1);


end
hold off;
AllTracks = [AllTracks, ConfirmedTracks];

figure;
hold on;
for k = 1:length(AllTracks)
    for l = 1:length(AllTracks(k).StateEstimateHistorySet)
        if l == 1
            estimated_x_position = zeros(1,length(AllTracks(k).StateEstimateHistorySet));
            estimated_y_position = zeros(1,length(AllTracks(k).StateEstimateHistorySet));
            estimated_x_velocity = zeros(1,length(AllTracks(k).StateEstimateHistorySet));
            estimated_y_velocity = zeros(1,length(AllTracks(k).StateEstimateHistorySet));
        end
    current_state = AllTracks(k).StateEstimateHistorySet{l};
    estimated_x_position(l) = current_state(1);
    estimated_y_position(l) = current_state(2);
    estimated_x_velocity(l) = current_state(3);
    estimated_y_velocity(l) = current_state(4);
    end
    plot(x_position_true,y_position_true,LineWidth=1.5,Color="#77AC30");
    plot(estimated_x_position,estimated_y_position,LineWidth=1.5,Color="#A2142F")
    title("True Target vs. Estimated Target Trajectory");
    ylabel("y position");
    ylim([0,yrange]);
    xlabel("x-position");
    xlim([0,xrange]);
    legend("True Target Trajectory","Estimated Target Trajectory");
    grid on;
end