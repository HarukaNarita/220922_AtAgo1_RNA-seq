library("stringr")
library(DESeq2)						#パッケージの読み込み

getDESeq2Result <- function(dir,control,treat) {
  # in_f <- paste0(dir, "/triplicate_" , control, "_", treat, ".txt")
  in_f <- paste0(dir, "/AtAgo1_RNA-seq_featureCounts_tRNA_" , control, "_", treat, ".txt")
  print(in_f)
  
  out_f1 <- paste0(str_sub(in_f,end = -5),"_DESeq2.txt")		#出力ファイル名を指定してout_f1に格納
  
  
  param_FDR <- 0.01					#false discovery rate (FDR)閾値を指定
  param_fig <- c(400, 380)				#ファイル出力時の横幅と縦幅を指定(単位はピクセル)
  
  count <- read.table(in_f, header=TRUE, row.names=1, sep="\t")	#in_fで指定したファイルの読み込み
  count <- as.matrix(count)
  
  dim(count)
  head(count)
  
  #前処理(DESeqDataSetオブジェクトの作成)
  # group <- data.frame(con = factor(c("A", "A", "A", "B", "B", "B")))
  group <- data.frame(con = factor(c("A", "A", "B", "B")))
  d <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~ con)	#DESeqDataSetオブジェクトdの作成
  
  #本番(DEG検出)
  d <- DESeq(d)									#DESeq2を実行
  tmp <- results(d)								#実行結果を抽出
  head(tmp)								 #中身確認
  baseMean <- tmp$baseMean             #baseMean
  log2FC <- tmp$log2FoldChange				 #log2FoldChangeをlog2FCに格納
  lfcSE <- tmp$lfcSE
  stat <- tmp$stat
  p.value <- tmp$pvalue								#p-valueをp.valueに格納
  p.value[is.na(p.value)] <- 1						#NAを1に置換している
  q.value <- tmp$padj								#adjusted p-valueをq.valueに格納
  q.value[is.na(q.value)] <- 1						#NAを1に置換している
  q.value_BH <- p.adjust(p.value, method="BH")
  ranking <- rank(p.value)							#p.valueでランキングした結果をrankingに格納
  sum(q.value < param_FDR)							#FDR閾値(q.value < param_FDR)を満たす遺伝子数を表示
  sum(q.value_BH < param_FDR)	#FDR閾値(q.value < param_FDR)を満たす遺伝子数を表示(TCCはp-valueをもとにBH法でq-value情報を得ている)
  
  #ファイルに保存(テキストファイル)
  tmp <- cbind(count, baseMean, log2FC, lfcSE, stat, p.value, q.value, q.value_BH, ranking)		#入力データの右側にDEG検出結果を結合したものをtmpに格納
  write.table(tmp, out_f1, sep="\t", append=F, quote=F, row.names=T)		#tmpの中身を指定したファイル名で保存
}
