
browns <- c("gold4", "goldenrod4", "orange4", "tan4")
cyans <- c("turquoise1", "cyan2", "turquoise2")
aquas <- c("cyan3","turquoise", "turquoise3", "mediumturquoise", "cadetblue3", "darkslategray3", "darkturquoise")
light_blues <- c("deepskyblue", "deepskyblue1", "deepskyblue2","steelblue1", "steelblue2", "skyblue", "skyblue1", "skyblue2", "skyblue3")
medium_blues <- c("steelblue", "steelblue3", "royalblue", "royalblue1", "royalblue2", "royalblue3", "cornflowerblue", "deepskyblue3", "dodgerblue", "dodgerblue1", "dodgerblue2", "dodgerblue3")
mints <- c("aquamarine", "aquamarine1", "aquamarine2")
turquoises <- c("cyan4", "lightseagreen", "turquoise4", "darkcyan")
dark_blues <- c("mediumblue", "darkblue", "blue", "blue1", "blue2", "blue3")
navy_blues <- c("royalblue4", "navy", "navyblue", "midnightblue", "darkslateblue", "blue4", "deepskyblue4", "dodgerblue4")
purples <- c("purple", "purple1", "purple2", "purple3", "blueviolet", "darkviolet", "darkorchid1", "darkorchid2", "darkorchid3", "darkorchid", "mediumorchid3")
plums <- c("plum", "plum1", "plum2", "plum3")
hot_pinks <- c("deeppink3", "violetred", "violetred1", "violetred2", "violetred3", "mediumvioletred", "deeppink", "deeppink1", "deeppink2", "maroon1", "maroon2", "maroon3")
dark_purples <-c ("purple4", "darkorchid4", "orchid4","slateblue4", "mediumorchid4", "mediumpurple4", "darkmagenta","magenta4")
orchids <- c("mediumorchid", "mediumorchid3", "mediumorchid1", "orchid", "orchid1", "orchid2", "orchid3")
magentas <- c("violet", "magenta", "magenta1", "magenta3", "magenta2")
lavenders <- c("slateblue3", "slateblue1", "slateblue2", "slateblue", "lightslateblue", "mediumpurple", "mediumpurple1", "mediumpurple2", "mediumpurple3")
bright_pinks <- c("palevioletred", "palevioletred1", "palevioletred2", "palevioletred3", "hotpink", "hotpink1", "hotpink2", "hotpink3")
pastel_pinks <- c("pink", "pink1", "pink2", "pink3", "lightpink1", "lightpink2", "lightpink3")
reds <- c("red", "red1", "red2", "red3", "brown", "brown1", "brown2", "firebrick", "firebrick1", "firebrick2", "firebrick3")
corals <- c("tomato", "tomato1", "tomato2", "tomato3", "salmon", "salmon1", "salmon2", "salmon3", "indianred3", "indianred2", "indianred1",  "coral", "coral1", "coral2", "lightcoral")
maroons <- c("violetred4", "red4", "palevioletred4", "deeppink4", "hotpink4", "firebrick4", "maroon4", "brown4", "darkred", "maroon")
golds <- c("orange3", "darkgoldenrod", "darkgoldenrod1", "darkgoldenrod2", "darkgoldenrod3", "gold3", "goldenrod", "goldenrod3")
greens <- c("darkgreen", "palegreen4", "darkseagreen4", "springgreen4", "chartreuse4", "forestgreen", "green4", "olivedrab", "olivedrab4", "darkolivegreen", "darkolivegreen4", "seagreen", "seagreen4")
soft_greens <- c("mediumaquamarine","aquamarine3", "springgreen3", "seagreen3", "palegreen3",  "mediumseagreen", "darkseagreen3")
yellows <- c("gold", "gold1", "gold2", "goldenrod1", "goldenrod2")
oranges <- c("tan1", "tan2", "orange", "orange1", "orange2", "chocolate1", "chocolate2", "darkorange", "darkorange1", "darkorange2")
red_oranges <- c("orangered", "orangered1", "orangered2", "orangered3")
bitter_greens <- c("green3", "limegreen", "olivedrab3", "yellowgreen", "chartreuse3")
light_greens <-c("lightgreen", "seagreen1", "seagreen2", "mediumspringgreen", "springgreen", "springgreen1", "springgreen2")
caramels <- c("chocolate", "chocolate3", "darkorange3")

colors <- list(dark_purples, hot_pinks, greens, dark_blues, yellows, maroons, browns, turquoises, orchids, pastel_pinks, corals, yellows, browns,  dark_blues, dark_purples, hot_pinks, greens, hot_pinks, dark_purples, dark_blues, yellows, lavenders, cyans, red_oranges, aquas, bitter_greens, navy_blues, browns, light_blues, maroons, mints, reds, oranges, bright_pinks, purples, medium_blues,caramels, plums, magentas, golds, soft_greens, light_greens)

plotColors <- function(colors){
  colors = sort(colors)
  numbers = seq(1, length(colors), 1)
  numbers2 = numbers
  df = data.frame(cbind(numbers, numbers2, colors))
  g = ggplot(df, aes(numbers, numbers2, color=colors))
  g = g + geom_point(size=5) + scale_color_manual(values=colors) + theme_classic()
  return(g)
}

generatePalette <- function(n, randomseed = 0){
  pal = vector(mode="character")
  for (i in 1:n){
    if (randomseed!=0){
      set.seed(randomseed)
      pal[i] = sample(colors[[i]], 1)
    } else {
      pal[i] = colors[[i]][1]
    }
  }
  return(pal)
}

