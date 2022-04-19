clc
clear
% Блок входных параметров для расчета
        simulationParams.horizontalElementsCount = 8;
        simulationParams.verticalElementsCount = 8;
        simulationParams.nUsers = 8;
        simulationParams.beamformerMethod = 'ZF';
        simulationParams.radAllocationMatrix = [];
% Запуск конструктора класса
beamformerObject = Beamformer(simulationParams);
% Расчет весовой матрицы
beamformerObject.getBeamformerWeights
% Скорость изменения фазы
phaseDistribution = angle(beamformerObject.beamformerWeights) * (180/pi);
