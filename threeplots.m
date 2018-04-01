function threeplots(Y1,Y2,Y3)
%CREATEFIGURE(Y1,Y2)
%  Y1:  vector of omega data
%  Y2:  vector of Tperc data

% Create figure
figure1 = figure;

% Create subplot
subplot1 = subplot(3,1,1,'Parent',figure1);
ylim(subplot1,[-.05 1.05]);
box(subplot1,'on');
hold(subplot1,'all');
% Create plot
plot(Y1,'Parent',subplot1);
% Create xlabel
xlabel('time');
% Create ylabel
ylabel('rational % of wealth');

% Create subplot
subplot2 = subplot(3,1,2,'Parent',figure1);
ylim(subplot2,[-.05 1.05]);
box(subplot2,'on');
hold(subplot2,'all');
% Create plot
plot(Y2,'Parent',subplot2);
% Create xlabel
xlabel('time');
% Create ylabel
ylabel('rational % of population');

% Create subplot
subplot3 = subplot(3,1,3,'Parent',figure1);
ylim(subplot3,[-.05 1.05]);
box(subplot3,'on');
hold(subplot3,'all');
% Create plot
plot(Y3,'Parent',subplot3);
% Create xlabel
xlabel('time');
% Create ylabel
ylabel('Gini coefficient');