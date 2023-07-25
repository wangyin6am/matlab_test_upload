clc
clear
%% 初始值
fs = 100;
space = fs * 0.4;
sleep_table = [0,1,2,3,5];
number_channel = 8;

% [sleep_flag, file_name] = xlsread('Z:\2_睡眠分析数据\乱七八糟的eeg格式\HRV分析进度\2-毕设脑电-1-72\前后半夜\毕设双向.xlsx');
[sleep_flag, file_name] = xlsread('Z:\2_睡眠分析数据\eeg格式\HRV分析进度\2-毕设脑电-1-72\前后半夜\毕设双向.xlsx');

% 1秒一个标签
% sleep_flag = sleep_flag(:,2);
file_name  = file_name(2:end,1);
% 滤波器
Hd = bandpass_window_100;
coeff_eegdata = Hd.Numerator;

%% 输出到文件夹
% Files = dir(fullfile('Z:\2_睡眠分析数据\乱七八糟的eeg格式\HRV分析进度\2-毕设脑电-1-72\数据'));
% FilesRR_all = ('Z:\2_睡眠分析数据\乱七八糟的eeg格式\HRV分析进度\2-毕设脑电-1-72\HRV结果\');
Files = dir(fullfile('Z:\2_睡眠分析数据\eeg格式\HRV分析进度\2-毕设脑电-1-72\数据'));
FilesRR_all = ('Z:\2_睡眠分析数据\eeg格式\HRV分析进度\2-毕设脑电-1-72\HRV结果\');
LengthFiles = length(Files);
%% a 
for j=3:LengthFiles
    RPoint_real_all= 0;
    name = Files(j).name;
    folder = Files(j).folder;
    % 加载睡眠分期
    FilesPath = [folder,'\',name,'\'];
    % 脑电数据
    FilesPathData = [folder,'\',name,'\','data_raw.mat'];
    FilesPathRR_all = [FilesRR_all,name,'_SleepHRV_StageOrder.xlsx'];
    FilesPathHypno = [FilesPath,'Hypnogram.mat'];
    FilesRRPath = [folder,'\',name,'\'];
    table_name = file_name{j-2};
    if(strcmp(table_name(2:end-1),name) == 0)
        disp([name,'名称不匹配']);
%         continue
    end
    % 保存RR间隔的文件
    FilesRR     = dir(strcat(FilesRRPath,'*.xlsx'));
    % RR间隔文件的数目
    FilesRR_number = length(FilesRR);
%     FilesPathOutFolder = [FilesPath,'输出结果\'];
%     mkdir(FilesPathOutFolder);
    % 读取结果文件
    RR_table_all  = readtable(FilesPathRR_all);
    RR_table_time = RR_table_all(:,1:4);
    % 按照阶段开始时间排序
    RR_table_time_sort = sortrows(RR_table_time,4);
    % 结果文件统计的阶段数目
    RR_table_number = numel(find(~isnan(RR_table_time_sort.Var4)));
    if(FilesRR_number ~= RR_table_number)
        disp([name,'RR文件表数目不匹配']);
        continue;
    end
    % 根据文件第一列统计阶段 W、N1、...
    RR_table_stage_name_cell = RR_table_time_sort(1:RR_table_number,1);
    RR_table_stage_name_all  = table2array(RR_table_stage_name_cell);
    flag = 0;
    % 循环读取RR间隔xlsx
    for jj = 1:FilesRR_number
        FileRR_one   = strcat('data_raw_RR-',num2str(jj),'.xlsx');
        % 不存在某个RR文件
        if(~exist([folder,'\',name,'\',FileRR_one],'file'))
            disp([name,FileRR_one,'RR间隔文件有问题']);
            flag = 1;
            break;
        end
        RR_table_one = readtable([folder,'\',name,'\',FileRR_one]);
        % 根据'_'分割
        RR_table_stage_name = strsplit(RR_table_stage_name_all{jj},'_');
        RR_table_stage_name = cell2mat(RR_table_stage_name(1));
        
        % W期
        if(strcmp(RR_table_stage_name,'W'))
            % 睡眠阶段
            stage.stage{jj}  = 0;
            % R峰时间
            stage.Rtime{jj}  = RR_table_one.RtimeR;
            % RR间隔
            stage.RR{jj}     = RR_table_one.RR;
            % R峰索引位置
            stage.Rpoint{jj} = round(RR_table_one.RtimeR * fs);
        end
        % N1期
        if(strcmp(RR_table_stage_name,'N1'))
            % 睡眠阶段
            stage.stage{jj}  = 1;
            % R峰时间
            stage.Rtime{jj}  = RR_table_one.RtimeR;
            % RR间隔
            stage.RR{jj}     = RR_table_one.RR;
            % R峰索引位置
            stage.Rpoint{jj} = round(RR_table_one.RtimeR * fs);
        end
        % N2期
        if(strcmp(RR_table_stage_name,'N2'))
            % 睡眠阶段
            stage.stage{jj}  = 2;
            % R峰时间
            stage.Rtime{jj}  = RR_table_one.RtimeR;
            % RR间隔
            stage.RR{jj}     = RR_table_one.RR;
            % R峰索引位置
            stage.Rpoint{jj} = round(RR_table_one.RtimeR * fs);
        end
        % N3期
        if(strcmp(RR_table_stage_name,'N3'))
            % 睡眠阶段
            stage.stage{jj}  = 3;
            % R峰时间
            stage.Rtime{jj}  = RR_table_one.RtimeR;
            % RR间隔
            stage.RR{jj}     = RR_table_one.RR;
            % R峰索引位置
            stage.Rpoint{jj} = round(RR_table_one.RtimeR * fs);
        end
        % REM期
        if(strcmp(RR_table_stage_name,'R'))
            % 睡眠阶段
            stage.stage{jj}  = 5;
            % R峰时间
            stage.Rtime{jj}  = RR_table_one.RtimeR;
            % RR间隔
            stage.RR{jj}     = RR_table_one.RR;
            % R峰索引位置
            stage.Rpoint{jj} = round(RR_table_one.RtimeR * fs);
        end
        
    end
    % RR间隔文件
    if(flag ==1)
        continue
    end

%%
    % 加载RR间隔文件，加载睡眠周期文件
    load(FilesPathHypno)
    eeg  = load(FilesPathData);
    data = eeg.eegdata;
    break_points = sleep_flag(j-2);
    hypno_second_before = hypno_second(1:break_points);
    hypno_second_after  = hypno_second(break_points+1:end);
    break_points = break_points * fs;
    %% 脑电减参考
    eeg_A1 = data(18,:);
    eeg_A2 = data(19,:);
    mean_A1A2 = (eeg_A1+eeg_A2) ./2;
    eeg = data(1:number_channel,:);
    eeg = eeg-mean_A1A2;
    clear mean_A1A2 eeg_A2 eeg_A1
    %% 滤波
    eegdata_filter_tmp = filtfilt(coeff_eegdata,1,eeg');
    eegdata_filter = eegdata_filter_tmp';
    clear eegdata_filter_tmp
    %% 找到五个阶段的时间段 + 伪迹
    %1.找到分期为5的地方，hypno_second(REM_idx_jump_idx(1,i):REM_idx_jump_idx(2,i)) 即为REM的时间片段，单位是秒
    Wake_idx = find(hypno_second == 0);
    Wake_idx = [Wake_idx,10];

    N1_idx = find(hypno_second == 1);
    N1_idx = [N1_idx,10];

    N2_idx = find(hypno_second == 2);
    N2_idx = [N2_idx,10];

    N3_idx = find(hypno_second == 3);
    N3_idx = [N3_idx,10];

    REM_idx = find(hypno_second == 5);
    REM_idx = [REM_idx,10];

    %2.找到发生跳变的地方
    Wake_idx_diff = diff(Wake_idx);
    N1_idx_diff = diff(N1_idx);
    N2_idx_diff = diff(N2_idx);
    N3_idx_diff = diff(N3_idx);
    REM_idx_diff = diff(REM_idx);
    tmp_Wake = [0 find(Wake_idx_diff ~= 1)];
    tmp_N1 = [0 find(N1_idx_diff ~= 1)];
    tmp_N2 = [0 find(N2_idx_diff ~= 1)];
    tmp_N3 = [0 find(N3_idx_diff ~= 1)];
    tmp_REM = [0 find(REM_idx_diff ~= 1)];
    clear *_diff
    % ———————
    tmp_Wake = tmp_Wake(ones(2,1),:);
    tmp_Wake(2,:) = tmp_Wake(2,:) + 1;
    tmp_N1 = tmp_N1(ones(2,1),:);
    tmp_N1(2,:) = tmp_N1(2,:) + 1;
    tmp_N2 = tmp_N2(ones(2,1),:);
    tmp_N2(2,:) = tmp_N2(2,:) + 1;
    tmp_N3 = tmp_N3(ones(2,1),:);
    tmp_N3(2,:) = tmp_N3(2,:) + 1;
    tmp_REM = tmp_REM(ones(2,1),:);
    tmp_REM(2,:) = tmp_REM(2,:) + 1;
    % ———————
    Wake_idx_jump_idx = zeros(2,size(tmp_Wake,2)-1);
    N1_idx_jump_idx = zeros(2,size(tmp_N1,2)-1);
    N2_idx_jump_idx = zeros(2,size(tmp_N2,2)-1);
    N3_idx_jump_idx = zeros(2,size(tmp_N3,2)-1);
    REM_idx_jump_idx = zeros(2,size(tmp_REM,2)-1);
    % ———————
    % hypno_second(REM_idx_jump_idx(1,i):REM_idx_jump_idx(2,i)) 即为某一周期的时间片段，单位是秒
    if(size(tmp_Wake,2) == 1)
        Wake_idx_jump_idx(1,1) = Wake_idx(1);
        Wake_idx_jump_idx(2,1) = Wake_idx(end);
    else
        for i=1:size(tmp_Wake,2)-1
            Wake_idx_jump_idx(1,i) = Wake_idx(tmp_Wake(2,i));
            Wake_idx_jump_idx(2,i) = Wake_idx(tmp_Wake(1,i+1));
        end
    end

    if(size(tmp_N1,2) == 1)
        N1_idx_jump_idx(1,1) = N1_idx(1);
        N1_idx_jump_idx(2,1) = N1_idx(end);
    else
        for i=1:size(tmp_N1,2)-1
            N1_idx_jump_idx(1,i) = N1_idx(tmp_N1(2,i));
            N1_idx_jump_idx(2,i) = N1_idx(tmp_N1(1,i+1));
        end
    end

    if(size(tmp_N2,2) == 1)
        N2_idx_jump_idx(1,1) = N2_idx(1);
        N2_idx_jump_idx(2,1) = N2_idx(end);
    else
        for i=1:size(tmp_N2,2)-1
            N2_idx_jump_idx(1,i) = N2_idx(tmp_N2(2,i));
            N2_idx_jump_idx(2,i) = N2_idx(tmp_N2(1,i+1));
        end
    end

    if(size(tmp_N3,2) == 1)
        N3_idx_jump_idx(1,1) = N3_idx(1);
        N3_idx_jump_idx(2,1) = N3_idx(end);
    else
        for i=1:size(tmp_N3,2)-1
            N3_idx_jump_idx(1,i) = N3_idx(tmp_N3(2,i));
            N3_idx_jump_idx(2,i) = N3_idx(tmp_N3(1,i+1));
        end
    end

    if(size(tmp_REM,2) == 1)
        REM_idx_jump_idx(1,1) = REM_idx(1);
        REM_idx_jump_idx(2,1) = REM_idx(end);
    else
        for i=1:size(tmp_REM,2)-1
            REM_idx_jump_idx(1,i) = REM_idx(tmp_REM(2,i));
            REM_idx_jump_idx(2,i) = REM_idx(tmp_REM(1,i+1));
        end
    end

    clear tmp*
    %% 匹配超过五分钟的数据
    % 超过五分钟的个数
    wn  = sum((Wake_idx_jump_idx(2,:) - Wake_idx_jump_idx(1,:)) >= 300);
    % 超过五分钟的索引
    fminWake_idx = find((Wake_idx_jump_idx(2,:) - Wake_idx_jump_idx(1,:)) >= 300);
    n1n = sum((N1_idx_jump_idx(2,:) - N1_idx_jump_idx(1,:)) >= 300);
    fminN1_idx = find((N1_idx_jump_idx(2,:) - N1_idx_jump_idx(1,:)) >= 300);
    n2n = sum((N2_idx_jump_idx(2,:) - N2_idx_jump_idx(1,:)) >= 300);
    fminN2_idx = find((N2_idx_jump_idx(2,:) - N2_idx_jump_idx(1,:)) >= 300);
    n3n = sum((N3_idx_jump_idx(2,:) - N3_idx_jump_idx(1,:)) >= 300);
    fminN3_idx = find((N3_idx_jump_idx(2,:) - N3_idx_jump_idx(1,:)) >= 300);
    remn = sum((REM_idx_jump_idx(2,:) - REM_idx_jump_idx(1,:)) >= 300);
    fminREM_idx = find((REM_idx_jump_idx(2,:) - REM_idx_jump_idx(1,:)) >= 300);
    % 数目匹配
    if((wn + n1n + n2n + n3n + remn) ~= RR_table_number)
        disp([name,'RR间隔5分钟阶段数目不匹配']);
        continue
    end
    % 根据5分钟片段分割脑电信号存入stage
    [wN,n1N,n2N,n3N,remN] = deal(0);
    for wk = 1:RR_table_number
        stage.HEP_sum{wk} = 0;
        switch stage.stage{wk}
            case 0
                wN = wN+1;
                n  = 0; 
                % 开始截止索引点
                stage.eegStartIdx{wk}   = Wake_idx_jump_idx(1,fminWake_idx(wN)) * fs;
                stage.eegEndIdx{wk}     = Wake_idx_jump_idx(2,fminWake_idx(wN)) * fs;
                % R峰出现的索引点
                stage.RPoint_real{wk}  = stage.Rpoint{wk} +  stage.eegStartIdx{wk} - 1;
                % R峰前后0.4秒的索引点
                stage.Rpoint_leftSpace{wk}  = stage.RPoint_real{wk} - 0.4 * fs;
                stage.Rpoint_rightSpace{wk} = stage.RPoint_real{wk} + 0.4 * fs;
            case 1
                n1N = n1N+1;
                n = 0;
                stage.eegStartIdx{wk}   = N1_idx_jump_idx(1,fminN1_idx(n1N)) * fs;
                stage.eegEndIdx{wk}     = N1_idx_jump_idx(2,fminN1_idx(n1N)) * fs;
                stage.RPoint_real{wk}  = stage.Rpoint{wk} +  stage.eegStartIdx{wk} - 1;
                stage.Rpoint_leftSpace{wk}  = stage.RPoint_real{wk} - 0.4 * fs;
                stage.Rpoint_rightSpace{wk} = stage.RPoint_real{wk} + 0.4 * fs;
            case 2
                n2N = n2N+1;
                n = 0;
                stage.eegStartIdx{wk}   = N2_idx_jump_idx(1,fminN2_idx(n2N)) * fs;
                stage.eegEndIdx{wk}     = N2_idx_jump_idx(2,fminN2_idx(n2N)) * fs;
                stage.RPoint_real{wk}  = stage.Rpoint{wk} +  stage.eegStartIdx{wk} - 1;
                stage.Rpoint_leftSpace{wk}  = stage.RPoint_real{wk} - 0.4 * fs;
                stage.Rpoint_rightSpace{wk} = stage.RPoint_real{wk} + 0.4 * fs;
            case 3
                n3N = n3N+1;
                n   = 0;
                stage.eegStartIdx{wk}   = N3_idx_jump_idx(1,fminN3_idx(n3N)) * fs;
                stage.eegEndIdx{wk}     = N3_idx_jump_idx(2,fminN3_idx(n3N)) * fs;
                stage.RPoint_real{wk}  = stage.Rpoint{wk} +  stage.eegStartIdx{wk} - 1;
                stage.Rpoint_leftSpace{wk}  = stage.RPoint_real{wk} - 0.4 * fs;
                stage.Rpoint_rightSpace{wk} = stage.RPoint_real{wk} + 0.4 * fs;
            case 5
                remN = remN+1;
                n = 0;
                stage.eegStartIdx{wk}   = REM_idx_jump_idx(1,fminREM_idx(remN)) * fs;
                stage.eegEndIdx{wk}     = REM_idx_jump_idx(2,fminREM_idx(remN)) * fs;
                stage.RPoint_real{wk}  = stage.Rpoint{wk} +  stage.eegStartIdx{wk} - 1;
                stage.Rpoint_leftSpace{wk}  = stage.RPoint_real{wk} - 0.4 * fs;
                stage.Rpoint_rightSpace{wk} = stage.RPoint_real{wk} + 0.4 * fs;
        % switch
        end
%         % 每个RR间隔文件中找到最合适的睡眠中间点
%         [~,mid_idx] = min(abs(break_points - stage.RPoint_real{wk}));
%         stage.mid_idx{wk}   = mid_idx;
%         stage.mid_point{wk} = stage.RPoint_real{wk}(mid_idx);

        % 整夜的R峰索引
        RPoint_real_all = [RPoint_real_all;stage.RPoint_real{wk}];

    end

    %% 前后半夜处理
    % 整夜R峰索引分割睡眠前后半夜
    [~,mid_idx] = min(abs(break_points - RPoint_real_all));
    mid_point   = RPoint_real_all(mid_idx);
    [stage.HEP.Wake_before_sum,stage.HEP.Wake_after_sum,stage.HEP.N1_before_sum,stage.HEP.N1_after_sum,stage.HEP.N2_before_sum,stage.HEP.N2_after_sum,stage.HEP.N3_before_sum,stage.HEP.N3_after_sum,stage.HEP.REM_before_sum,stage.HEP.REM_after_sum] = deal(0);
        % 根据5分钟片段分割脑电信号存入stage
    [wNL,n1NL,n2NL,n3NL,remNL,wNR,n1NR,n2NR,n3NR,remNR] = deal(0);
    [wN,n1N,n2N,n3N,remN] = deal(0);
    for wk = 1:RR_table_number
        [HEP_Wake_before,HEP_Wake_after,HEP_N1_before,HEP_N1_after,HEP_N2_before,HEP_N2_after,HEP_N3_before,HEP_N3_after,HEP_REM_before,HEP_REM_after] = deal(0);
        
        switch stage.stage{wk}
            case 0
                stage.HEP_Wake_before{wk} = 0;
                stage.HEP_Wake_after{wk}  = 0;
                wN = wN+1;
                nb  = 0; 
                na  = 0;
                for wkk = 1:size(stage.Rpoint_leftSpace{wk},1)
                    % 前半夜
                    if(stage.RPoint_real{wk}(wkk) < mid_point)
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        nb = nb + 1;
                        stage.before_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_Wake_before      = HEP_Wake_before + stage.before_eegdata{wk}{wkk};
                        end
                    % 后半夜
                    else
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        na = na + 1;
                        stage.after_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_Wake_after      = HEP_Wake_after + stage.after_eegdata{wk}{wkk};
                        end
                    end
                end
                num.Wake.na{wN} = na;
                num.Wake.nb{wN} = nb;
                stage.HEP.Wake_before_sum = stage.HEP.Wake_before_sum + HEP_Wake_before;
                stage.HEP.Wake_after_sum  = stage.HEP.Wake_after_sum + HEP_Wake_after;
            case 1
                n1N = n1N+1;
                stage.HEP_N1_before{wk} = 0;
                stage.HEP_N1_after{wk}  = 0;
                nb  = 0; 
                na  = 0;
                for wkk = 1:size(stage.Rpoint_leftSpace{wk},1)
                    % 前半夜
                    if(stage.RPoint_real{wk}(wkk) < mid_point)
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        nb = nb + 1;
                        stage.before_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_N1_before       = HEP_N1_before + stage.before_eegdata{wk}{wkk};
                        end
                    % 后半夜
                    else
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        na = na + 1;
                        stage.after_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_N1_after      = HEP_N1_after + stage.after_eegdata{wk}{wkk};
                        end
                    end
                    
                end
                num.N1.na{n1N} = na;
                num.N1.nb{n1N} = nb;
                stage.HEP.N1_before_sum = stage.HEP.N1_before_sum + HEP_N1_before;
                stage.HEP.N1_after_sum  = stage.HEP.N1_after_sum + HEP_N1_after;
            case 2
                stage.HEP_N2_before{wk} = 0;
                stage.HEP_N2_after{wk}  = 0;
                n2N = n2N+1;
                nb  = 0; 
                na  = 0;
                for wkk = 1:size(stage.Rpoint_leftSpace{wk},1)
                    % 前半夜
                    if(stage.RPoint_real{wk}(wkk) < mid_point)
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        nb = nb + 1;
                        stage.before_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_N2_before = HEP_N2_before + stage.before_eegdata{wk}{wkk};
                        end
                    % 后半夜
                    else
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        na = na + 1;
                        stage.after_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_N2_after       = HEP_N2_after + stage.after_eegdata{wk}{wkk};
                        end
                    end
                    
                end
                num.N2.na{n2N} = na;
                num.N2.nb{n2N} = nb;
                stage.HEP.N2_before_sum = stage.HEP.N2_before_sum + HEP_N2_before;
                stage.HEP.N2_after_sum  = stage.HEP.N2_after_sum + HEP_N2_after;
            case 3
                stage.HEP_N3_before{wk} = 0;
                stage.HEP_N3_after{wk}  = 0;
                n3N = n3N+1;
                nb  = 0; 
                na  = 0;
                for wkk = 1:size(stage.Rpoint_leftSpace{wk},1)
                    % 前半夜
                    if(stage.RPoint_real{wk}(wkk) < mid_point)
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        nb = nb + 1;
                        stage.before_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_N3_before = HEP_N3_before + stage.before_eegdata{wk}{wkk};
                        end
                    % 后半夜
                    else
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        na = na + 1;
                        stage.after_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_N3_after       = HEP_N3_after + stage.after_eegdata{wk}{wkk};
                        end
                    end
                    
                end
                num.N3.na{n3N} = na;
                num.N3.nb{n3N} = nb;
                stage.HEP.N3_before_sum = stage.HEP.N3_before_sum + HEP_N3_before;
                stage.HEP.N3_after_sum  = stage.HEP.N3_after_sum + HEP_N3_after;
            case 5
                stage.HEP_REM_before{wk} = 0;
                stage.HEP_REM_after{wk}  = 0;
                remN = remN+1;
                nb  = 0; 
                na  = 0;
                for wkk = 1:size(stage.Rpoint_leftSpace{wk},1)
                    % 前半夜
                    if(stage.RPoint_real{wk}(wkk) < mid_point)
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        nb = nb + 1;
                        stage.before_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_REM_before       = HEP_REM_before + stage.before_eegdata{wk}{wkk};
                        end
                    % 后半夜
                    else
                        if(stage.Rpoint_leftSpace{wk}(wkk)>=1 && stage.Rpoint_rightSpace{wk}(wkk)<=size(eegdata_filter,2))
                        na = na + 1;
                        stage.after_eegdata{wk}{wkk}  = eegdata_filter(:,stage.Rpoint_leftSpace{wk}(wkk):stage.Rpoint_rightSpace{wk}(wkk));
                        HEP_REM_after       = HEP_REM_after + stage.after_eegdata{wk}{wkk};
                        end
                    end
                    
                end
                num.REM.na{remN} = na;
                num.REM.nb{remN} = nb;
                stage.HEP.REM_before_sum = stage.HEP.REM_before_sum + HEP_REM_before;
                stage.HEP.REM_after_sum  = stage.HEP.REM_after_sum + HEP_REM_after;
        % switch
        end
%         % 每个RR间隔文件中找到最合适的睡眠中间点
%         [~,mid_idx] = min(abs(break_points - stage.RPoint_real{wk}));
%         stage.mid_idx{wk}   = mid_idx;
%         stage.mid_point{wk} = stage.RPoint_real{wk}(mid_idx);

    end
    if(wN ~= 0)
        % HEP=总/n
        stage.HEP.Wake_before = stage.HEP.Wake_before_sum / sum(cell2mat(num.Wake.nb));
        stage.HEP.Wake_after  = stage.HEP.Wake_after_sum / sum(cell2mat(num.Wake.na));
    % 如果没有，置零
    else
        stage.HEP.Wake_before = zeros(number_channel,0.8 * fs + 1);
        stage.HEP.Wake_after  = zeros(number_channel,0.8 * fs + 1);
    end
    if(n1N ~= 0)
        % HEP=总/n
        stage.HEP.N1_before = stage.HEP.N1_before_sum / sum(cell2mat(num.N1.nb));
        stage.HEP.N1_after  = stage.HEP.N1_after_sum / sum(cell2mat(num.N1.na));
    else
        stage.HEP.N1_before = zeros(number_channel,0.8 * fs + 1);
        stage.HEP.N1_after  = zeros(number_channel,0.8 * fs + 1);
    end
    if(n2N ~= 0)
        % HEP=总/n
        stage.HEP.N2_before = stage.HEP.N2_before_sum / sum(cell2mat(num.N2.nb));
        stage.HEP.N2_after  = stage.HEP.N2_after_sum / sum(cell2mat(num.N2.na));
    else
        stage.HEP.N2_before = zeros(number_channel,0.8 * fs + 1);
        stage.HEP.N2_after  = zeros(number_channel,0.8 * fs + 1);
    end
    if(n3N ~= 0)
        % HEP=总/n
        stage.HEP.N3_before = stage.HEP.N3_before_sum / sum(cell2mat(num.N3.nb));
        stage.HEP.N3_after  = stage.HEP.N3_after_sum / sum(cell2mat(num.N3.na));
    else
        stage.HEP.N3_before = zeros(number_channel,0.8 * fs + 1);
        stage.HEP.N3_after  = zeros(number_channel,0.8 * fs + 1);
    end
    if(remN ~= 0)
        % HEP=总/n
        stage.HEP.REM_before = stage.HEP.REM_before_sum / sum(cell2mat(num.REM.nb));
        stage.HEP.REM_after  = stage.HEP.REM_after_sum / sum(cell2mat(num.REM.na));
    else
        stage.HEP.REM_before = zeros(number_channel,0.8 * fs + 1);
        stage.HEP.REM_after  = zeros(number_channel,0.8 * fs + 1);
    end

    %%
    
    people.stage{j-2,1}  = stage;
    people.name{j-2,1} = name;
    clear stage
    disp(j-2);


%wy
end
%% 将所有数据存储为特定矩阵
for j=3:LengthFiles
%     HEP.Wake_before(:,:,j-2) = people.stage{j-2,1}.HEP.Wake_before;
%     HEP.Wake_after(:,:,j-2)  = people.stage{j-2,1}.HEP.Wake_after;
%     HEP.N1_before(:,:,j-2)   = people.stage{j-2,1}.HEP.N1_before;
%     HEP.N1_after(:,:,j-2)    = people.stage{j-2,1}.HEP.N1_after;
    HEP.N2_before(:,:,j-2)   = people.stage{j-2,1}.HEP.N2_before;
    HEP.N2_after(:,:,j-2)    = people.stage{j-2,1}.HEP.N2_after;
%     HEP.N3_before(:,:,j-2)   = people.stage{j-2,1}.HEP.N3_before;
%     HEP.N3_after(:,:,j-2)    = people.stage{j-2,1}.HEP.N3_after;
%     HEP.REM_before(:,:,j-2)  = people.stage{j-2,1}.HEP.REM_before;
%     HEP.REM_after(:,:,j-2)   = people.stage{j-2,1}.HEP.REM_after;
end
% HEP.Wake_all = cat(3,HEP.Wake_before,HEP.Wake_after);
% HEP.N1_all = cat(3,HEP.N1_before,HEP.N1_after);
HEP.N2_all = cat(3,HEP.N2_before,HEP.N2_after);
% HEP.N3_all = cat(3,HEP.N3_before,HEP.N3_after);
% HEP.REM_all = cat(3,HEP.REM_before,HEP.REM_after);

%% 将HH:MM:SS时间转换为秒数时间
function results = calcTime(time)
        str = strsplit(time,':');
        x1 = corr(str2double(str(1)));
        x2 = corr(str2double(str(2)));
        x3 = corr(str2double(str(3)));
        results = x1*3600+x2*60+x3;
end
%% 带通滤波器
function Hd = bandpass_window_100
    %BANDPASS_WINDOW_100 Returns a discrete-time filter object.
    
    %
    % MATLAB Code
    % Generated by MATLAB(R) 7.14 and the Signal Processing Toolbox 6.17.
    %
    %
    
    % FIR Window Bandpass filter designed using the FIR1 function.
    
    % All frequency values are in Hz.
    Fs = 100;  % Sampling Frequency
    
    N    = 300;     % Order
    Fc1  = 0.3;      % First Cutoff Frequency
    Fc2  = 35;       % Second Cutoff Frequency
    flag = 'scale';  % Sampling Flag
    Beta = 0.5;       % Window Parameter
    % Create the window vector for the design algorithm.
    win = kaiser(N+1, Beta);
    % Calculate the coefficients using the FIR1 function.
    b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);  
    Hd = dfilt.dffir(b);
end





