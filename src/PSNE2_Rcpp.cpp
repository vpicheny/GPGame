// #include <Rcpp.h>
// using namespace Rcpp;
//
//
//
//
// // Computes the pure strategy Nash Equilibria in a 2-player finite game.
// // ---------------------------------------------------------------------
// // input:
// //
// // NS(i) : total number of pure strategies of player (i)
// //
// // Poff : list of two matrices of dimension (NS[1],NS[2])
// // Poff(i,is1,is2) : payoff matrix of player (i) when
// //                   player (1) plays is1 and player (2) plays is2
// //
// // output : NE
// //          NE(1:2,j) is the Nash equilibrium number j
// //
// // ---------------------------------------------------------------------
// // A. Habbal Inria 2016-05
// // V. Picheny INRA 2016-06
// // M. Binois Chicago Booth 2016-06
// // [[Rcpp::export]]
// IntegerVector PSNE_Rcpp(NumericVector NS, NumericMatrix Poff1, NumericMatrix Poff2){
//   IntegerMatrix iB1(NS(0), NS(1));
//   IntegerMatrix iB2(NS(0), NS(1));
//   std::vector<int> NE;
//   double B1, B2;
//
//   // Set of best replies of player 2 to moves by player 1
//   for(int is1 = 0; is1 < NS(0); is1++ ){
//     B2 = min(Poff2(is1,_));
//     for(int is2 = 0; is2 < NS(1); is2++){
//       if(Poff2(is1, is2) == B2){
//         iB1(is1, is2) = 1;
//       }
//     }
//   }
//
//   // for (is1 in 1:NS[1]){
//   //   B2 <- min(Poff[[2]][is1,1:NS[2]])
//   //   for (is2 in 1:NS[2]) {
//   //     if( Poff[[2]][is1,is2]==B2 ){
//   //       iB1[is1,is2] <- 1
//   //     }
//   //   }
//   // }
//
//   // Set of best replies of player 1 to moves by player 2
//   for(int is2 = 0; is2 < NS(1); is2++){
//     B1 = min(Poff1(_, is2));
//     for(int is1 = 0; is1 < NS(0); is1++){
//       if(Poff1(is1, is2) == B1){
//         iB2(is1, is2) = 1;
//       }
//     }
//   }
//
//   // for (is2 in 1:NS[2]){
//   //   B1 <- min(Poff[[1]][1:NS[1],is2])
//   //   for (is1 in 1:NS[1]) {
//   //     if(Poff[[1]][is1,is2]==B1) {
//   //       iB2[is1,is2] <- 1
//   //     }
//   //   }
//   // }
//
//   // Enumerate pure strategy Nash equilibria
//   // NE <- c()
//   for(int is1 = 0; is1 < NS(0); is1++){
//     for(int is2 = 0; is2 < NS(1); is2++){
//       if(iB1(is1, is2) + iB2(is1, is2) == 2){
//         NE.push_back(is1+1);
//         NE.push_back(is2+1);
//       }
//     }
//   }
//
//   // for (is1 in 1:NS[1]){
//   //   for (is2 in 1:NS[2]) {
//   //     if (iB1[is1,is2]+iB2[is1,is2] == 2) {
//   //       NE <- rbind(NE, c(is1, is2))
//   //     }
//   //   }
//   // }
//
//   return(as<IntegerVector>(wrap(NE)));
//
// }
//
//
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
//
// /*** R
// library(microbenchmark)
//
// ## To have the same entries and outputs than PSNE/PSNE2
// PSNE3 <- function(NS, Poff){
//   tmp <- PSNE_Rcpp(NS, Poff[[1]], Poff[[2]])
//   if(length(tmp) == 0)
//     return(NULL)
//   return(matrix(tmp, length(tmp)/2, byrow = T))
// }
//
// Z <- matrix(response.grid, ncol=nobj)
// Z1 <- matrix(Z[,1], n.grid, n.grid)
// Z2 <- matrix(Z[,2], n.grid, n.grid)
// Poff <- list(Z1, Z2)
//
//
// microbenchmark(PSNE(NS=rep(n.grid, 2), Poff = Poff),
//                PSNE2(NS=rep(n.grid, 2), Poff = Poff),
//                PSNE3(rep(n.grid, 2), Poff = Poff))
//
// */
