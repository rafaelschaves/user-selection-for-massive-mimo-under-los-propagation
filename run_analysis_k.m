MC = 10000;                                                                % Size of the outer Monte Carlo ensemble (Varies the channel realizarions)

snr_db = 10;

channel_type = 'ur-los';

for M = [64 128 256]                                                       % Number of antennas at the base station
    for K = [1 2 4 8 16 32 64 128 256 512 1024]                            % Number of users at the cell
        run analysis_k
    end
end