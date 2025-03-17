classdef WaveformAnalyzer < handle
    %% Описание класса
    %
    % 1. Класс читает данные (во временной области) на выходе OFDM модулятора сигнала, а также информацию о параметрах формирователя
    %
    % 2. Строит метрики: спектральная плотность мощности в частотной области, графическое представление созвездия на комплексной плоскости,
    % среднеквадратичное значение модуля вектора ошибки (EVM)
    %
    % Входные данные:
    %
    % waveformSource - массив содержащий отчеты baseband сигнала во временной области на выходе OFDM модулятора
    %
    % waveformInfo - структура с параметрами OFDM модулятора и пейлоуда:
    %       Nfft               - кол-во спектрально-временных отчетов дискретного преобразования Фурье
    %       SampleRate         - частота семплирования [Гц]
    %       CyclicPrefixLengths/SymbolLengths - длины циклического преффикса и OFDM символов [кол-во временных отчетов]
    %       SymbolsCount       - кол-во символов на слот радиокадра
    %       subCarriersCount   - кол-во поднесущих
    %       payloadSymbols     - информационные символы
    %       payloadSymbolsIdxs - индексы ресурсных элементов отведенные для передачи payloadSymbols
    %
    % Поля класса:
    %
    %       rmsEvm            - среднеквадратичное значение модуля вектора ошибки
    %       waveformMeanPower - среднеквадратичное значение мощности сигнала
    %       channelBandwidth  - ширина полосы канала
    %       noiseMeanPower    - среднеквадратичное значение мощности шума
    %       modulationType    - тип модуляционной схемы
    %       waveformDuration  - длина анализируемого сигнала
    %
    properties (Access = private)
        nfft
        sampleRate
        cyclicPrefixLengths
        symbolLengths
        windowing
        symbolPhases
        symbolsPerSlot
        symbolsCount
        payloadSymbols
        subCarriersCount
        payloadSymbolsdxs
        rxWaveform
        sampleCount
        powerSpectrumDensity
        f
    end


    properties
        rmsEvm
        waveformMeanPower
        channelBandwidth
        noiseMeanPower
        modulationType
        waveformDuration
        dopplerShift
        modulationOrder
        
    end

    methods
        function this = WaveformAnalyzer(waveformInfo,waveformSource)
            % Конструктор класса. Чтение waveform-ы во временной области и структуры с информацией
            % необходимой для дальнейшей обработки данных и заполнения полей класса
            this.nfft = waveformInfo.Nfft;
            this.sampleRate = waveformInfo.SampleRate;
            this.cyclicPrefixLengths = waveformInfo.CyclicPrefixLengths;
            this.symbolLengths = waveformInfo.SymbolLengths;
            this.windowing = waveformInfo.Windowing;
            this.symbolPhases = waveformInfo.SymbolPhases;
            this.symbolsPerSlot = waveformInfo.SymbolsPerSlot;
            this.symbolsCount = waveformInfo.symbolsCount;
            this.payloadSymbols = waveformInfo.payloadSymbols;
            this.subCarriersCount = waveformInfo.subCarriersCount;
            this.payloadSymbolsdxs = waveformInfo.payloadSymbolsIdxs;
            this.rxWaveform = waveformSource;
        end

        function calcWaveformParameters(this)
            %Метод класса реализующий расчет параметров сигнала
            %Результат:
            %среднеквадратичное значение мощности сигнала
            %длина анализируемого сигнала
            %позиционность модуляции
            %ширина полосы канала

            %Среднеквадратичное значение мощности сигнала
            this.waveformMeanPower = rms(this.rxWaveform);

            %Число очтетов сигнала
            this.sampleCount = length(this.rxWaveform);
            %Расчет длительности сигнала
            this.waveformDuration = length(this.rxWaveform)/this.sampleRate;

            %Определение уникальных значений payload
            uniquePayloadSymbols = unique (this.payloadSymbols);
            %Определение позиционности модуляции
            this.modulationOrder = length(uniquePayloadSymbols);

            %Расчет СПМ
            this.powerSpectrumDensity = fftshift(1/(this.sampleRate*this.sampleCount)*abs(fft(this.rxWaveform)).^2);
            %Диапаизон значений СПМ
            this.f = linspace (-this.sampleRate/2,this.sampleRate/2,this.sampleCount);
            %Среднее значение СПМ
            averagePowerSpectrumDensity = mean(this.powerSpectrumDensity);

            %Определение начала и конца полосы сигнала
            indBandwidthFirst = find (this.powerSpectrumDensity>averagePowerSpectrumDensity,1,'first');
            indBandwidthLast = find (this.powerSpectrumDensity>averagePowerSpectrumDensity,1,'last');
            %Расчет ширины полосы сигнала
            this.channelBandwidth = this.f(indBandwidthLast)-this.f(indBandwidthFirst);

        end

        function calcDopplerShift(this)
            %Метод класса реализующий расчет доплеровского смещения
            %частоты Результат:
            %Вектор смещения частоты символов OFDM сигнала

            ofdmSlot214 = zeros(this.symbolLengths(2),this.symbolsCount-1);

            %Цикл определения OFDM символов
            for i = 1:this.symbolsCount
                %первый символ
                if i==1
                    ofdmSlot1 = this.rxWaveform(1:this.symbolLengths(1));
                %Символы со 2 по 14
                else
                    ofdmSlot214 (:,i-1) =  this.rxWaveform(this.symbolLengths(1)+1+this.symbolLengths(2)*(i-2):...
                        this.symbolLengths(1)+this.symbolLengths(2)*(i-1));
                end
            end

            %БПФ символов
            dataSubcar1 = fft(ofdmSlot1);
            dataSubcar214 = fft(ofdmSlot214);

            %Значения циклического префикса и окончания сигнала 1ого
            %символа
            dataCyclicPrefix1 = dataSubcar1(1:this.cyclicPrefixLengths(1));
            dataEndSlot1 = dataSubcar1(end-this.cyclicPrefixLengths(1)+1);

            %Значения циклического префикса и окончания сигнала для
            %символов со 2 по 14
            dataCyclicPrefix214 = zeros(this.cyclicPrefixLengths(2),this.symbolsCount-1);
            dataEndSlot214 = dataCyclicPrefix214;
            for i=1:this.symbolsCount-1
                dataCyclicPrefix214(:,i) = dataSubcar214(1:this.cyclicPrefixLengths(2),i);
                dataEndSlot214(:,i) = dataSubcar214(end-this.cyclicPrefixLengths(2)+1,i);
            end
            
            %Цикл определения доплеровского смещения частоты каждого OFDM
            %символа
            for i=1:this.symbolsCount
                %Первый символ
                if i==1
                    %корреляционная ф-ция циклического префикса и конца
                    %символа
                    [corrFunc,~] = xcorr(dataEndSlot1,dataCyclicPrefix1);
                    %индекс максимума корреляционной ф-ции
                    [~,ind_max] = max(abs(corrFunc));
                    %Определение частоты смещения
                    dopplerPhase = angle(corrFunc(ind_max));
                    this.dopplerShift(i) = dopplerPhase/(2*pi*(this.nfft/this.sampleRate));
                %символы со 2 по 14
                else
                    %корреляционная ф-ция циклического префикса и конца
                    %символа
                    [corrFunc,~] = xcorr(dataEndSlot214(:,i-1),dataCyclicPrefix214(:,i-1));
                    %индекс максимума корреляционной ф-ции
                    [~,ind_max] = max(abs(corrFunc));
                    %Определение частоты смещения
                    dopplerPhase = angle(corrFunc(ind_max));
                    this.dopplerShift(i) = dopplerPhase/(2*pi*(this.nfft/this.sampleRate));
                end
            end
        end

        function plotPowerSpectrumDensity(this)
            %Метод реализующий графическое построение спектральной
            %плотности мощности
            figure
            plot(this.f,pow2db(this.powerSpectrumDensity));
            title('Cпектральная плотность мощности');
            xlabel('Частота, Гц');
            ylabel('Мощность/Частоту, дБ/Гц');
            grid on;
        end

        function plotPayloadConstellation(this)
            %метод реализующий графическое построение сигнального созвездия
            %Payload
            scatterplot (this.payloadSymbols);
            title('Payload');
            grid on;
        end

        function calcEvmPerformance(this)

        end
    end
end