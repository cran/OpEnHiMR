pred_hmc <- function(fasta_file_path){
  fasta_file_test <- fasta_file_path
  ################################Training##############################
  # Function to pad DNA sequence
  pad_dna_sequence <- function(sequence, desired_length, padding_char = "N") {
    current_length <- nchar(sequence)

    if (current_length >= desired_length) {
      # If the sequence is equal or longer than the desired length, no padding needed
      return(sequence)
    } else {
      # Calculate the number of characters to pad
      pad_length <- desired_length - current_length

      # Pad the sequence with the specified character
      padded_sequence <- paste0(sequence, strrep(padding_char, pad_length))

      return(padded_sequence)
    }
  }

  # Function to pad DNA sequences in a multifasta file
  pad_multifasta_file <- function(input_file, desired_length, padding_char = "N") {
    # Read multifasta file
    fasta_sequences <- readDNAStringSet(input_file, format = "fasta")

    # Pad each sequence
    padded_sequences <- sapply(as.character(fasta_sequences), function(seq) {
      pad_dna_sequence(seq, desired_length, padding_char)
    })

    # Convert the padded sequences to a data frame
    padded_sequences_df <- data.frame(ID = names(fasta_sequences), Sequence = padded_sequences, stringsAsFactors = FALSE)

    return(padded_sequences_df)
  }

  #######################Fasta to Tabular format##################################

  FastaToTabular <- function (filename){

    #read fasta file

    file1 <- readLines(filename)

    #find the genename location by grepping >

    location <- which((str_sub(file1,1,1))==">")


    name=c()
    sequence =c()
    for ( i in 1:length(location)){
      name_line = location[i]
      name1 = file1[name_line]
      name=c(name,name1)
      start= location[i]+1
      end = location[i+1]-1
      if ( i < length (location)){

        end=end

      } else {

        end=length(file1)
      }

      lines = start:end
      sequence1= as.character(paste(file1[lines],collapse = ""))
      sequence =c(sequence,sequence1)
    }

    #now create table using name and sequence vector

    data <- tibble(name,sequence)



    #finally export the file
    #before that remove preexisting file
    unlink(c("dna_table.csv"),force=TRUE)
    as.matrix(data,"dna_table.csv")

    #function ends
  }
  #########################alphabetcheck###########################
  alphabetCheck<-function (sequences, alphabet = "aa", label = c())
  {
    if (length(sequences) == 0) {
      stop("ERROR: sequence parameter is empty")
    }
    if (length(label) != 0 && length(label) != length(sequences)) {
      stop("ERROR: The lenght of the label vector and the number of sequences do not match!")
    }
    if (alphabet == "rna") {
      alphabet <- c("A", "C", "G", "U", "N")
    }
    else if (alphabet == "dna") {
      alphabet <- c("A", "C", "G", "T", "N")
    }
    else if (alphabet == "aa") {
      alphabet <- c("A", "C", "D", "E", "F", "G", "H", "I",
                    "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V",
                    "W", "Y")
    }
    else {
      stop("ERROR: alphabet shoud be 'dna' or 'rna' or 'aa' ")
    }
    alphabetCheck = sapply(sequences, function(i) all(strsplit(i,
                                                               split = "")[[1]] %in% alphabet))
    flag = 0
    if (length(label) == length(sequences)) {
      flag = 1
      label = label[alphabetCheck]
    }
    else if (length(label) > 0 && length(label) != length(sequences)) {
      stop("ERROR: The number of labels is not equal to the number of sequences!")
    }
    if (is.null(names(sequences))) {
      names(sequences) <- as.character(1:length(sequences))
    }
    nonstanSeq <- names(sequences)[!alphabetCheck]
    if (length(nonstanSeq) != 0) {
      nonstanSeq <- toString(nonstanSeq)
      warMessage <- paste("The sequences (", nonstanSeq, ") were deleted. They contained non-standard alphabets")
      message(warMessage)
    }
    sequences = sequences[alphabetCheck]
    if (length(sequences) == 0) {
      stop("All sequences contained non-standard alphabets. No sequences remained for analysis :) ")
    }
    if (flag == 1) {
      names(label) = names(sequences)
    }
    seq_lab <- list(sequences = sequences, Lab = label)
    return(seq_lab)
  }
  #################################NCP_DNA############################

  ncp_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                     label = c())
  {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- fa.read(seqs, alphabet = "dna")
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else if (is.vector(seqs)) {
      seqs <- sapply(seqs, toupper)
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else {
      stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
    }
    lenSeqs <- sapply(seqs, nchar)
    nucs <- list(A = c(1, 1, 1), C = c(0, 0, 1), G = c(1, 0,
                                                       0), T = c(0, 1, 0), U = c(0, 1, 0), N = c(0,0,0))
    numSeqs <- length(seqs)
    if (outFormat == "mat") {
      if (length(unique(lenSeqs)) > 1) {
        stop("ERROR: All sequences should have the same length in 'mat' mode.
           For sequences with different lengths, please use 'txt' for outFormat parameter")
      }
      if (binaryType == "strBin") {
        nucs <- c(A = "111", C = "001", G = "100", T = "010",
                  U = "010", N = "000")
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        colnames(featureMatrix) <- paste("ncp_pos", 1:lenSeqs[1],
                                         sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "logicBin") {
        nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                    TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                   FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("ncp_pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "numBin") {
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("ncp_pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else {
        stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
      }
      return(featureMatrix)
    }
    else if (outFormat == "txt") {
      nucs <- c(A = "111", C = "001", G = "100", T = "010",
                U = "010", N = "000")
      counter <- 0
      namesSeqs <- names(seqs)
      codes <- lapply(seqs, function(x) {
        counter <- counter + 1
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        namecods <- namesSeqs[counter]
        cods <- unlist(cods)
        cods <- c(namecods, cods)
        temp <- paste(cods, collapse = "\t")
        write(temp, outputFileDist, append = TRUE)
      })
    }
    else {
      stop("ERROR: outFormat should be 'mat' or 'txt' ")
    }
  }

  #######MBE_DNA############


  mbe_dna<-function (seqs, binaryType = "numBin", outFormat = "mat", outputFileDist = "",
                     label = c())
  {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- fa.read(seqs, alphabet = "dna")
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else if (is.vector(seqs)) {
      seqs <- sapply(seqs, toupper)
      seqs_Lab <- alphabetCheck(seqs, alphabet = "dna", label)
      seqs <- seqs_Lab[[1]]
      label <- seqs_Lab[[2]]
    }
    else {
      stop("ERROR: Input sequence is not in the correct format. It should be a FASTA file or a string vector.")
    }
    lenSeqs <- sapply(seqs, nchar)
    nucs <- list(A = c(1, 0, 0, 0), C = c(0, 0, 1, 0), G = c(0, 0, 0,
                                                             1), T = c(0, 1, 0, 0), U = c(0, 1, 0, 0), N = c(0,0,0,0))
    numSeqs <- length(seqs)
    if (outFormat == "mat") {
      if (length(unique(lenSeqs)) > 1) {
        stop("ERROR: All sequences should have the same length in 'mat' mode.
           For sequences with different lengths, please use 'txt' for outFormat parameter")
      }
      if (binaryType == "strBin") {
        nucs <- c(A = "1000", C = "0010", G = "0001", T = "0100",
                  U = "0100", N = "0000")
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        colnames(featureMatrix) <- paste("mbe_pos", 1:lenSeqs[1],
                                         sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "logicBin") {
        nucs <- list(A = c(TRUE, TRUE, TRUE), C = c(FALSE,
                                                    TRUE, FALSE), G = c(TRUE, FALSE, FALSE), T = c(FALSE,
                                                                                                   FALSE, TRUE), U = c(FALSE, FALSE, TRUE))
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("mbe_pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else if (binaryType == "numBin") {
        featureMatrix <- sapply(seqs, function(x) {
          charList <- unlist(strsplit(x, split = ""))
          cods <- nucs[charList]
          cods <- unlist(cods)
          return(cods)
        })
        featureMatrix <- t(featureMatrix)
        temp1 <- rep(c("P", "A", "H"), lenSeqs[1])
        temp2 <- rep(1:lenSeqs[1], each = 3)
        colnames(featureMatrix) <- paste("mbe_pos", temp2, "-",
                                         temp1, sep = "")
        row.names(featureMatrix) <- names(seqs)
      }
      else {
        stop("ERROR! Choose one of 'strBin', 'logicBin', or 'numBin' for binaryFormat")
      }
      return(featureMatrix)
    }
    else if (outFormat == "txt") {
      nucs <- c(A = "1000", C = "0010", G = "0001", T = "0100",
                U = "0100", N = "0000")
      counter <- 0
      namesSeqs <- names(seqs)
      codes <- lapply(seqs, function(x) {
        counter <- counter + 1
        charList <- unlist(strsplit(x, split = ""))
        cods <- nucs[charList]
        namecods <- namesSeqs[counter]
        cods <- unlist(cods)
        cods <- c(namecods, cods)
        temp <- paste(cods, collapse = "\t")
        write(temp, outputFileDist, append = TRUE)
      })
    }
    else {
      stop("ERROR: outFormat should be 'mat' or 'txt' ")
    }
  }

  ###########################GC_Content##########################
  GC.content <- function(fasta_file){
    x <- read.fasta(file=fasta_file)
    tt<-function(x){
      res<-GC(x)
      val=round(res,4)
      return(val)
    }

    f_res<-lapply(x,tt)
    s=data.frame(f_res)

    rownames(s) <- c("GC-content")

    w=t(s)
    return(w)
  }
  ################################ONF#################################
  oligo.freq <- function(fasta_file,f){
    x<- readDNAStringSet(fasta_file)
    y <- oligonucleotideFrequency(x,width = f)
    z <- data.frame(y)
    rownames(z) <- names(x)

    return(z)
  }

  ##################model calling########################
  svm_model <- readRDS(url("https://github.com/snehaiasri/OpEnHiMR/raw/main/model/svm.rds"))
  rf_model <- readRDS(url("https://github.com/snehaiasri/OpEnHiMR/raw/main/model/rf.rds"))
  gb_model <- readRDS(url("https://github.com/snehaiasri/OpEnHiMR/raw/main/model/gb.rds"))
  ###########################Testing#################################
  ## Sequence Padding
  seq_object_test <- read.fasta(fasta_file_test)
  lengths_test <- getLength(seq_object_test)
  ids_test <- names(seq_object_test)
  result_table_test <- data.frame(Sequence_ID = ids_test, Length = lengths_test)
  desired_length_test <- max(result_table_test$Length) + 400
  padded_sequences_test <- pad_multifasta_file(fasta_file_test, desired_length_test)

  ##Feature Extraction
  data_test <- as.vector(padded_sequences_test[,2])
  ncp_mat_test <- as.matrix(ncp_dna(seqs = data_test,binaryType="strBin",outFormat="mat"))
  sequence_test <- rownames(ncp_mat_test)
  seq_id_test <- padded_sequences_test[,1]
  ncp_test <- cbind(seq_id_test,sequence_test,ncp_mat_test)
  rownames(ncp_test) <- seq_id_test
  ncp_temp_test <- data.frame(ncp_test[,-1], stringsAsFactors = FALSE)
  ncp_final_test <- as.data.frame(apply(ncp_temp_test[,-1], 2, as.numeric))
  mbe_test <- as.matrix(mbe_dna(seqs = data_test,binaryType="strBin",outFormat="mat"))
  mbe_temp_test <- data.frame(mbe_test, stringsAsFactors = FALSE)
  mbe_final_test <- as.data.frame(apply(mbe_temp_test, 2, as.numeric))
  log_gc_temp_test<-log((GC.content(fasta_file_test))*100, base = exp(1))
  log_gc_test<-as.data.frame(as.numeric((ifelse(log_gc_temp_test>0,log_gc_temp_test,'0'))))
  onf2_test <- oligo.freq(fasta_file_test, 2)
  onf3_test <- oligo.freq(fasta_file_test, 3)
  temp1_test <- cbind(onf2_test, onf3_test, gcc =log_gc_test[,1], ncp_final_test, mbe_final_test)
  inputData_test <-as.data.frame(temp1_test)
  zero_columns <- colSums(inputData_test) == 0
  inputData_test <- inputData_test[, !zero_columns]

  #############################Feature_Selection##########################

  colname_rf_test <- c("TAT", "ncp_pos1", "GC", "TA", "mbe_pos1", "TTA", "GGC", "TAA", "ATT", "AA", "gcc",
                       "TGC", "AAA", "mbe_pos2", "TGG", "TAC", "AT", "AGC", "AAT", "TT", "GCT", "CA", "TTT",
                       "CTG", "CC", "TAG", "ATA", "GCA", "GCC", "GG", "CGC", "mbe_pos385", "AGT", "ACG",
                       "mbe_pos130", "CGT", "CGG", "GCG", "AAC", "mbe_pos352", "CG", "CGA", "mbe_pos237",
                       "ACT", "CAG","CTC", "mbe_pos149")
  test_data <- inputData_test[, colname_rf_test]


  ############################Prediction########################
  predicted_prob_svm_temp <- predict(svm_model, newdata = test_data, type = "prob", probability = TRUE)
  predicted_prob_svm <- attr(predicted_prob_svm_temp, "probabilities")
  predicted_value_svm <- predict(svm_model, newdata = test_data)
  predicted_prob_rf <- predict(rf_model, newdata = test_data, type = "prob")
  predicted_value_rf <- predict(rf_model, newdata = test_data)
  predicted_prob_gb <- predict(gb_model, newdata = test_data, type = "prob")
  predicted_value_gb <- predict(gb_model, newdata = test_data)
  test_res_en_prob <-  cbind(SVM = predicted_prob_svm, RF = predicted_prob_rf, GB = predicted_prob_gb)


  ##################Ensemble###################
  # Define weights
  weights <- c(SVM = 0.0936, RF = 0.8129, GB = 0.0935)

  # Apply weights
  weighted_data <- as.data.frame(test_res_en_prob)
  weighted_data[1:4] <- test_res_en_prob[1:4] * weights["SVM"]
  weighted_data[5:8] <- test_res_en_prob[5:8] * weights["RF"]
  weighted_data[9:12] <- test_res_en_prob[9:12] * weights["GB"]

  weighted_data$Class_1 <- rowSums(weighted_data[c("SVM.1", "RF.1", "GB.1")])
  weighted_data$Class_2 <- rowSums(weighted_data[c("SVM.2", "RF.2", "GB.2")])
  weighted_data$Class_3 <- rowSums(weighted_data[c("SVM.3", "RF.3", "GB.3")])
  weighted_data$Class_4 <- rowSums(weighted_data[c("SVM.4", "RF.4", "GB.4")])


  weighted_data$Predicted_Class <- apply(weighted_data[, c("Class_1", "Class_2", "Class_3", "Class_4")], 1, which.max)


  final_results <- weighted_data[, c("Class_1", "Class_2", "Class_3", "Class_4", "Predicted_Class")]

  final_results$Predicted_Class_Label <- ifelse(final_results$Predicted_Class == 1, "H3K27me3",
                                                ifelse(final_results$Predicted_Class == 2, "H3K9ac",
                                                       ifelse(final_results$Predicted_Class == 3, "H3K4me3",
                                                              ifelse(final_results$Predicted_Class == 4, "No Modification", NA))))
  final_pred <- data.frame(Ids = rownames(final_results),
                           Modification = final_results$Predicted_Class_Label,
                           Probability = round(apply(final_results[, c("Class_1", "Class_2", "Class_3", "Class_4")], 1, max), 2))
  rownames(final_pred) <- NULL
  return(final_pred)
}
