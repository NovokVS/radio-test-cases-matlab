% Скрипт для запуска waveformAnalyxer
clc
clear
addpath waveform/

%Импорт данных
        waveformInfo = importdata('waveformInfo.mat');
        waveformSource = importdata('waveformSource.mat');
%Запуск конструктора класса
waveformAnalyzerObject = WaveformAnalyzer(waveformInfo,waveformSource);
%Расчет параметров сигнала
waveformAnalyzerObject.calcWaveformParameters;
%Построение сигнального созвездие Payload
waveformAnalyzerObject.plotPayloadConstellation;
%Построение графика спектральной плотности мощности
waveformAnalyzerObject.plotPowerSpectrumDensity;
%Расчет доплеровского смещения частоты
waveformAnalyzerObject.calcDopplerShift;