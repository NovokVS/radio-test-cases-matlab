classdef Beamformer < handle
    %% Описание класса
    % 1. Класс производит генерацию канальных коэффициентов системы MIMO
    % порядка [nTransmitAntennas x nUsers].
    %
    % 2. На основе канальных коэффициентов в частотной области производит
    % расчет матрицы прекодирования используя методы "MRT" или "ZF".
    %
    % 3. Согласно произведенным вычислениям рассчитывается ресурсная
    % производительность радиопередачи - спектральная эффективность.
    %
    % 4. КА несет на борту многоэлементную АФАР с количеством элементов
    % по горизонтали "horizontalElementsCount" и по вертикали
    % "verticalElementsCount".
    %
    % 5. Абонентское устройство имеет однолучевую антенну на основе
    % параболической зеркальной антенны.

    properties (Access = private)
        horizontalElementsCount
        verticalElementsCount
        allocationMatrix
        multiuserLayout
        nUsers
        snrdB
    end

    properties
        beamformerWeights
        beamformerMethod
        channelCoeffs
        spectralPerformance
    end

    methods

        function this = Beamformer(simulationParams)
            % Конструктор класса
            % Парсинг входной структуры
            this.horizontalElementsCount = simulationParams.horizontalElementsCount;
            this.verticalElementsCount = simulationParams.verticalElementsCount;
            this.nUsers = simulationParams.nUsers;
            this.beamformerMethod = simulationParams.beamformerMethod;
            this.allocationMatrix = simulationParams.radAllocationMatrix;

            % Генерация канальных коэффициентов
            this.calcChannelRealization;

            % Расчет матриц прекодирования
            this.calcBeamformerWeights;

            % Расчет спектральной эффективности радиопередачи с учетом
            % матрицы прекодирования
            this.calcSpectralPerformance;
        end

        function calcChannelRealization(this)
            % Метод класса, реализующий получение канальных коэффициентов в
            % частотной области. Результат:
            % матрица [nUsers x (horizontalElementsCount x verticalElementsCount)]

            % Симуляционный объект
            s = qd_simulation_parameters;
            % Частота на которой генерируются канальные коэффициенты [11.7 ГГц]
            s.center_frequency = 11.7e9;
            % Выключаем строку выполнения
            s.show_progress_bars = 0;
            % Объект распположения
            layout = qd_layout(s);
            % Расположим антенну на КА на высоте 500 км
            layout.tx_position = [0 0 500e3]';
            % Повернем антенну КА таким образом чтобы она смотрела в надир
            layout.tx_track.orientation = [0; -pi/2; 0];
            % Зададим количество пользователей
            layout.no_rx = this.nUsers;
            % Произвольно расположим пользователей на Земле в пределе 100 км от
            % подспутниковой точки КА
            layout.randomize_rx_positions(100e3, 0, 0, 0);
            % Определеим ориаентации антенн пользовательских терминалов
            uePos = layout.tx_track.initial_position;
            orientation = zeros(3, 1);
            for userIdx = 1 : this.nUsers
                uePosTemp = layout.rx_track(1, userIdx).initial_position;
                rt = uePos - uePosTemp;
                rt = rt / norm(rt);
                orientation(2) = asin(rt(3));
                orientation(3) = atan2(rt(2), rt(1));
                layout.rx_track(1, userIdx).orientation = orientation;
            end
            % Задаем антенную систему КА (многоэлементная решетка)
            % Поляризационный режим
            polInd = 1;
            % Межэлементное расстояние (пол длины волны)
            elementSpacing = 0.5;
            % Количество панелей
            verticalPanels = 1;
            horizontalPanels = 1;
            % Межпанельное расстояние
            vertPanSpacing = 0;
            horizPanSpacing = 0;
            % Передающая антенна КА
            satTransmitAnt = qd_arrayant('3gpp-nr', this.horizontalElementsCount,...
                this.verticalElementsCount,...
                s.center_frequency, polInd, 0, elementSpacing, verticalPanels,...
                horizontalPanels, vertPanSpacing, horizPanSpacing);
            % Приемные антенны пользователей
            ueReceiveAnt = qd_arrayant('parabolic', 0.3,  s.center_frequency,...
                [], 1, 1, [], []);
            % Назначаем антенны пользователям и КА
            layout.tx_array = satTransmitAnt;
            layout.rx_array = ueReceiveAnt;
            % Определеим сценарий
            layout.set_scenario('Freespace');
            % Генерация канальных коэффициентов
            channelBuilder = layout.init_builder(1, 1, 'accurate');
            channelBuilder.gen_parameters;
            c = channelBuilder.get_channels;
            % Сырые канальные коэффициенты
            rawChannelCoeffs = transpose(...
                reshape([ c.coeff ], [ satTransmitAnt.no_elements, layout.no_rx ]));
            % Среденее значение коэффициента передачи
            averacalcransmissionGain = sum(sum(abs(rawChannelCoeffs) .^ 2))...
                / (satTransmitAnt.no_elements * layout.no_rx);
            % Нормированная матрица канала связи
            this.channelCoeffs = rawChannelCoeffs ./ sqrt(averacalcransmissionGain);
            this.multiuserLayout = layout;
        end % Конец function calcChannelRealization

        function calcBeamformerWeights(this)
            % Метод класса, реализующий расчет весовых коэффициентов для
            % преобразования передаваемой информации. Можно использовать два метода
            % расчета "MRT" или "ZF". Результат:
            % матрица [(horizontalElementsCount x verticalElementsCount) x nUsers]

            % Количество передающих антенн
            nTransmitAntennas = size(this.channelCoeffs, 2);
            % Если все передающие антенны могут работать со всеми пользователями то
            if isempty(this.allocationMatrix)
                D = repmat(eye(nTransmitAntennas), [ 1 1 this.nUsers ]);
            end
            % Предопределим матрицу под весовые коэффициенты
            this.beamformerWeights = zeros(size(this.channelCoeffs'));
            switch this.beamformerMethod
                case 'MRT'
                    % Расссчитаем весовые коэффициенты матрицы
                    for userIdx = 1 : this.nUsers
                        % коэффициенты канала c учетом тензора D
                        channelVector = (this.channelCoeffs(userIdx, :) * D(:, :, userIdx))';
                        % нормируем результат
                        this.beamformerWeights(:, userIdx) = channelVector / norm(channelVector);
                    end
                case 'ZF'
                    % Расссчитаем весовые коэффициенты матрицы
                    for userIdx = 1 : this.nUsers
                        % коэффициенты канала c учетом тензора D
                        effectiveChannel = (this.channelCoeffs * D(:, :, userIdx))';
                        % Расчет псевдоинверсии канальной матрицы
                        channelInversion = effectiveChannel / (effectiveChannel' * effectiveChannel);
                        % Нормирование весовых коэффициентов
                        this.beamformerWeights(:,userIdx) = channelInversion(:, userIdx) / norm(channelInversion(:, userIdx));
                    end
                otherwise
                    error('Выбранный тип матрицы прекодирования не найден!');
            end
        end % Конец function calcBeamformerWeights(this)

        function calcSpectralPerformance(this)
            % Метод класса, реализующий расчет спектральной эффективности
            % радиопередачи с учетом матрицы прекодирования полученной в расчете
            % и канальной матрицы системы MIMO

            % Задание диапазона рассматриваемых значений ОСШ в dB
            this.snrdB = linspace(- 40,5,45);
            % Расчет принятого сигнала
            channelGains = abs(this.channelCoeffs * this.beamformerWeights) .^ 2;
            % Расчет полезной составляющей сигнала (прямого канала прохождения)
            signalGains = diag(channelGains);
            % Расчет интерференционной составляющей сигнала
            % (сумма всех элементов строки матрицы за вычетом диагонального элемента)
            interferenceGains = sum(channelGains, 2) - signalGains;
            % Расчет спектральной эффективности по формуле Шеннона для k-го пользователя 
            % и для системы (суммированием по всем строкам/пользователям)
            this.spectralPerformance = sum(log2(1 + signalGains ./ (db2pow(- this.snrdB) + interferenceGains)), 1);
        end % Конец function calcSpectralPerformance

        function vuzailizeLayout(this)
            % Метод класса, реализующий графическое изображение расчетного случая.
            % Расположения КА и пользователей будут представлены в декартовой
            % системе координат.
            this.multiuserLayout.visualize();
        end

        function vuzailizeSpectralPerformance(this)
            % Метод класса, реализующий графическое изображение расчетного случая.
            % Выводятся зависимости спектральной эффективности от ОСШ
            legg = {0};
            for objIdx=1:numel(this)
                plot(this(objIdx).snrdB, this(objIdx).spectralPerformance);
                legg(objIdx) = {this(objIdx).beamformerMethod};
                hold on
            end
            hold off
            legend(legg);
            title('Precoding scheme comparision');
            xlabel('SNR, dB');
            ylabel('Spectral Performance, bits/s/Hz');
            xticks('auto');
            yticks('auto');
            grid on;
        end

    end
end

