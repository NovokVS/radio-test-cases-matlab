clc
clear
addpath quadriga_src/

% Блок 1 входных параметров для расчета
        simulationParams1.horizontalElementsCount = 8;
        simulationParams1.verticalElementsCount = 8;
        simulationParams1.nUsers = 8;
        simulationParams1.beamformerMethod = 'MRT';
        simulationParams1.radAllocationMatrix = [];
% Запуск конструктора класса 1
beamformerObject1 = Beamformer(simulationParams1);
% Расчет канальных коэффициентов
beamformerObject1.calcChannelRealization
% Расчет матриц прекодирования
beamformerObject1.calcBeamformerWeights
% Расчет спектральной эффективности
beamformerObject1.calcSpectralPerformance

% Блок 2 входных параметров для расчета
        simulationParams2.horizontalElementsCount = 8;
        simulationParams2.verticalElementsCount = 8;
        simulationParams2.nUsers = 8;
        simulationParams2.beamformerMethod = 'ZF';
        simulationParams2.radAllocationMatrix = [];
% Запуск конструктора класса 2
beamformerObject2 = Beamformer(simulationParams2);
% Расчет канальных коэффициентов
beamformerObject2.calcChannelRealization
% Расчет матриц прекодирования
beamformerObject2.calcBeamformerWeights
% Расчет спектральной эффективности
beamformerObject2.calcSpectralPerformance

% Создание массива объектов из разных блоков входных параметров
beamformerObjects = [beamformerObject1, beamformerObject2];

% Вывод зависимостей спектральной эффективности от ОСШ
beamformerObjects.vuzailizeSpectralPerformance
