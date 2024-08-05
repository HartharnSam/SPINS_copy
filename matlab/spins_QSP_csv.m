function [T1_max, T1_min, S1_max, S1_min, pdfCount] = spins_QSP_csv(filename)
  %% Reads in the csv file and interprets it in the correct way for the user.
  %%
  %% Input:
  %%    Filename (string): name of csv file.
  %%
  %% Output:
  %%    T1_max (float): upper limit of the maximum-valued bin of T1 tracer
  %%    T1_min (float): lower limit of the minimum-valued bin of T1 tracer
  %%    S1_max (float): upper limit of the maximum-valued bin of S1 tracer
  %%    S1_min (float): lower limit of the minimum-valued bin of S1 tracer
  A = importdata(filename);
  T1_max = A(1, 1);
  T1_min = A(1, 2);
  S1_max = A(1, 3);
  S1_min = A(1, 4);
  pdfCount = A(2:end, :);
  pdfCount = pdfCount / sum(sum(pdfCount));
end
