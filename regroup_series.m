function e = regroup_series()

e_2p5mm = load('e_2p5mm.mat');
e_2p0mm = load('e_2p0mm.mat');

e = e_2p5mm.e.copyObject();
e.serie = (e.serie(:) + e_2p0mm.e.serie(:));

e.reorderSeries('name');

end