

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0002376","immune system process", 0.600,-1.426,-0.678, 4.886,-1.5720,0.986,0.000),
c("GO:0006950","response to stress", 4.575, 6.190, 3.436, 5.769,-26.4782,0.687,0.000),
c("GO:0019748","secondary metabolic process", 0.138,-0.100, 8.194, 4.247,-9.8440,0.902,0.000),
c("GO:0050896","response to stimulus",12.210,-0.765,-4.221, 6.195,-21.7788,0.987,0.000),
c("GO:0051704","multi-organism process", 0.751,-2.478,-3.902, 4.984,-10.7510,0.986,0.000),
c("GO:0009813","flavonoid biosynthetic process", 0.016,-3.435,-1.130, 3.317,-2.5440,0.923,0.009),
c("GO:0042445","hormone metabolic process", 0.090, 0.142,-5.428, 4.064,-2.6054,0.918,0.010),
c("GO:0006790","sulfur compound metabolic process", 1.822,-4.101,-3.375, 5.369,-6.1013,0.936,0.013),
c("GO:0009812","flavonoid metabolic process", 0.018,-2.058,-2.799, 3.355,-2.0060,0.960,0.026),
c("GO:0071554","cell wall organization or biogenesis", 0.950, 2.941,-5.375, 5.086,-1.4800,0.959,0.038),
c("GO:0007568","aging", 0.088, 2.575,-4.646, 4.052,-1.8406,0.944,0.060),
c("GO:0044262","cellular carbohydrate metabolic process", 1.257, 0.591,-2.919, 5.208,-2.2421,0.934,0.071),
c("GO:0019742","pentacyclic triterpenoid metabolic process", 0.000,-4.290,-2.105, 0.845,-1.5878,0.892,0.074),
c("GO:0005975","carbohydrate metabolic process", 5.260,-0.587,-1.766, 5.829,-4.6593,0.942,0.078),
c("GO:0019760","glucosinolate metabolic process", 0.003,-3.878, 5.305, 2.625,-8.2123,0.700,0.093),
c("GO:0042430","indole-containing compound metabolic process", 0.285,-3.178, 2.736, 4.563,-6.6784,0.867,0.175),
c("GO:0055114","oxidation-reduction process",15.060,-1.589, 7.939, 6.286,-3.6800,0.854,0.185),
c("GO:0052544","defense response by callose deposition in cell wall", 0.001, 6.094,-0.132, 1.944,-2.4514,0.720,0.246),
c("GO:1901605","alpha-amino acid metabolic process", 3.625,-2.843, 6.690, 5.668,-6.7588,0.696,0.321),
c("GO:0044106","cellular amine metabolic process", 0.446,-3.692, 2.801, 4.757,-4.7968,0.796,0.324),
c("GO:0009308","amine metabolic process", 0.521,-4.673, 2.669, 4.825,-4.4579,0.875,0.330),
c("GO:0033037","polysaccharide localization", 0.047, 4.820,-4.815, 3.784,-1.5878,0.960,0.347),
c("GO:0042180","cellular ketone metabolic process", 0.423,-1.623, 8.504, 4.735,-6.1013,0.822,0.365),
c("GO:0016137","glycoside metabolic process", 0.031,-3.590, 4.682, 3.603,-5.3178,0.857,0.371),
c("GO:0009611","response to wounding", 0.127, 6.396, 1.935, 4.212,-12.7788,0.726,0.379),
c("GO:0016143","S-glycoside metabolic process", 0.003,-4.762, 4.098, 2.625,-8.2123,0.847,0.381),
c("GO:0051707","response to other organism", 0.299, 6.646, 3.110, 4.584,-15.2148,0.664,0.415),
c("GO:0016138","glycoside biosynthetic process", 0.024,-4.825, 4.672, 3.490,-2.1398,0.813,0.416),
c("GO:0009607","response to biotic stimulus", 0.342, 5.107, 3.864, 4.643,-15.6820,0.744,0.421),
c("GO:0006970","response to osmotic stress", 0.082, 6.003, 2.006, 4.022,-3.9650,0.688,0.427),
c("GO:1901700","response to oxygen-containing compound", 0.503, 6.743, 4.178, 4.810,-16.1252,0.657,0.441),
c("GO:0009719","response to endogenous stimulus", 0.526, 6.452, 3.846, 4.829,-9.8599,0.736,0.443),
c("GO:1901657","glycosyl compound metabolic process", 2.961,-2.668, 6.119, 5.580,-1.3366,0.825,0.446),
c("GO:0009628","response to abiotic stimulus", 0.571, 5.268, 4.274, 4.865,-7.7206,0.734,0.448),
c("GO:0009744","response to sucrose", 0.007, 5.732, 5.677, 2.952,-2.0391,0.698,0.454),
c("GO:0046417","chorismate metabolic process", 0.240,-1.946, 8.501, 4.489,-2.1381,0.793,0.473),
c("GO:0006722","triterpenoid metabolic process", 0.005,-4.850, 0.004, 2.775,-1.3223,0.883,0.483),
c("GO:0010817","regulation of hormone levels", 0.161, 0.473,-5.484, 4.314,-2.1398,0.951,0.486),
c("GO:0006952","defense response", 0.568, 6.366, 2.444, 4.863,-11.7026,0.697,0.497),
c("GO:0009605","response to external stimulus", 1.370, 5.671, 3.524, 5.245,-17.0951,0.716,0.501),
c("GO:0010193","response to ozone", 0.002, 7.455, 2.651, 2.427,-1.5878,0.694,0.517),
c("GO:0044272","sulfur compound biosynthetic process", 1.235,-4.629, 1.279, 5.200,-2.9136,0.844,0.534),
c("GO:0072330","monocarboxylic acid biosynthetic process", 0.940,-3.397, 6.780, 5.081,-1.6002,0.728,0.554),
c("GO:0043648","dicarboxylic acid metabolic process", 1.019,-2.343, 7.985, 5.116,-2.9250,0.768,0.560),
c("GO:0006082","organic acid metabolic process", 9.086,-2.415, 7.628, 6.067,-8.5148,0.763,0.561),
c("GO:0042221","response to chemical", 3.071, 5.913, 3.408, 5.595,-18.8648,0.697,0.562),
c("GO:0001101","response to acid chemical", 0.124, 6.378, 4.974, 4.201,-13.8837,0.687,0.563),
c("GO:0046677","response to antibiotic", 0.128, 6.165, 4.883, 4.215,-2.2946,0.687,0.565),
c("GO:0006979","response to oxidative stress", 0.575, 6.604, 2.362, 4.868,-9.1913,0.696,0.570),
c("GO:1901698","response to nitrogen compound", 0.178, 6.565, 4.762, 4.359,-1.3223,0.680,0.581),
c("GO:0006855","drug transmembrane transport", 0.189, 6.923, 2.234, 4.384,-1.7122,0.632,0.584),
c("GO:0042493","response to drug", 0.266, 6.655, 4.552, 4.534,-7.3650,0.672,0.602),
c("GO:0080167","response to karrikin", 0.006, 4.477, 2.758, 2.874,-1.8983,0.753,0.609),
c("GO:0010035","response to inorganic substance", 0.317, 6.032, 4.508, 4.609,-2.9321,0.668,0.611),
c("GO:0009725","response to hormone", 0.335, 6.640, 4.415, 4.633,-9.8809,0.641,0.614),
c("GO:0009073","aromatic amino acid family biosynthetic process", 0.540,-3.538, 5.982, 4.841,-1.5819,0.704,0.638),
c("GO:1901606","alpha-amino acid catabolic process", 0.569,-2.931, 7.052, 4.863,-2.9321,0.702,0.642),
c("GO:0009625","response to insect", 0.001, 4.464, 1.581, 1.991,-2.6054,0.758,0.658),
c("GO:0052386","cell wall thickening", 0.001, 2.594,-3.455, 2.117,-2.0066,0.942,0.659),
c("GO:0009072","aromatic amino acid family metabolic process", 0.719,-2.943, 6.948, 4.965,-3.3182,0.737,0.660),
c("GO:0071456","cellular response to hypoxia", 0.022, 6.988, 2.985, 3.457,-3.2236,0.628,0.662),
c("GO:0010200","response to chitin", 0.004, 5.543, 5.692, 2.755,-2.6541,0.705,0.664),
c("GO:0052545","callose localization", 0.002, 4.564,-4.924, 2.318,-1.6424,0.947,0.667),
c("GO:0002213","defense response to insect", 0.001, 6.186, 0.805, 2.104,-2.0912,0.728,0.668),
c("GO:0019745","pentacyclic triterpenoid biosynthetic process", 0.000,-5.068,-0.724, 0.602,-1.5878,0.869,0.669),
c("GO:0010033","response to organic substance", 0.900, 6.228, 4.190, 5.062,-12.5211,0.642,0.675),
c("GO:0033554","cellular response to stress", 2.967, 6.583, 2.822, 5.581,-1.7204,0.633,0.679),
c("GO:0031667","response to nutrient levels", 0.151, 5.168, 2.976, 4.288,-3.4853,0.688,0.688),
c("GO:0009651","response to salt stress", 0.043, 6.143, 1.669, 3.741,-3.3549,0.700,0.691),
c("GO:0009743","response to carbohydrate", 0.041, 6.147, 5.313, 3.726,-4.2383,0.671,0.692),
c("GO:0009851","auxin biosynthetic process", 0.001,-0.119,-5.354, 2.185,-2.5320,0.880,0.700));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
