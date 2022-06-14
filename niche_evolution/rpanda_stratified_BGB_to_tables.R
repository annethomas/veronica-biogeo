#getAnywhere(.stratified_BGB_to_tables)
# A single object matching '.stratified_BGB_to_tables' was found
# It was found in the following places
# namespace:RPANDA
# with value


stratified_BGB_to_tables <- function (tree, clado_events_tables, i) 
{
  hold <- events_txt_list_into_events_table(events_txt_list = clado_events_tables[[i]]$anagenetic_events_txt_below_node, 
                                             trtable = NULL, recalc_abs_ages = TRUE)
  rows <- which(clado_events_tables[[i]][, 34] != "none")
  if (length(rows) < 1) {
    ana.int <- NULL
  }
  else {
    row.len <- length(rows)
    if (row.len < dim(hold)[1]) {
      ho <- clado_events_tables[[i]][, 34]
      torep <- vector()
      for (n in 1:length(ho)) {
        torep <- c(torep, length(strsplit(ho[n], ";")[[1]]))
      }
      torep <- torep[rows]
      rows2 <- vector()
      for (m in 1:row.len) {
        rows2 <- c(rows2, rep(rows[m], torep[m]))
      }
      rows <- rows2
    }
    desc <- clado_events_tables[[i]][rows, 1]
    anc <- clado_events_tables[[i]][rows, 5]
    anmat <- cbind(desc, anc, hold)
    m <- regexpr("->", anmat$event_txt)
    regs <- regmatches(as.character(anmat$event_txt), m, 
                       invert = T)
    new_rangetxt <- unlist(lapply(regs, function(x) x[2]))
    nodenum_at_top_of_branch <- anmat[, 1]
    abs_event_time <- anmat$abs_event_time
    ana.int <- data.frame(abs_event_time, nodenum_at_top_of_branch, 
                          new_rangetxt)
    ana.int$new_rangetxt <- as.character(ana.int$new_rangetxt)
  }
  clado.int <- clado_events_tables[[i]]
  clado.int$clado_event_txt <- as.character(clado.int$clado_event_txt)
  if (any(is.na(clado.int$clado_event_txt))) {
    clado.int$clado_event_txt[which(is.na(clado.int$clado_event_txt))] <- ""
  }
  clado.int <- clado.int[which(clado.int$clado_event_txt != 
                                 ""), ]
  clado.int[, 12] <- clado.int$label
  node_ht = round(max(nodeHeights(tree)) - clado.int$time_bp, 
                  5)
  clado.int[, 9] <- node_ht
  lef <- vector()
  rig <- vector()
  for (i in 1:length(clado.int$node)) {
    nin <- clado.int$node[i]
    hol <- tree$edge[tree$edge[, 1] == nin, 2]
    lef <- c(lef, hol[1])
    rig <- c(rig, hol[2])
  }
  clado.int[, 15] <- lef
  clado.int[, 16] <- rig
  clado.int[, 20] <- clado.int$clado_event_txt
  colnames(clado.int)[c(9, 12, 15, 16, 20)] <- c("node_ht", 
                                                 "label", "left_desc_nodes", "right_desc_nodes", 
                                                 "clado_event_txt")
  int2 <- matrix(nrow = length(tree$tip.label), ncol = dim(clado.int)[2])
  int2 <- as.data.frame(int2)
  int2[, 1] <- 1:length(tree$tip.label)
  int2[, 12] <- tree$tip.label
  colnames(int2) <- colnames(clado.int)
  int <- rbind(clado.int, int2)
  clado.int <- int
  return(list(ana.int = ana.int, clado.int = clado.int))
}

events_txt_list_into_events_table<-function (events_txt_list, trtable = NULL, recalc_abs_ages = TRUE) 
{
  if (is.null(events_txt_list)) {
    errortxt = paste("\nWARNING in events_txt_list_into_events_table(): your events_txt_list has NO events!\n\nThis means your tree has NO d/e/a events across the whole tree.\nThis is *expected* e.g. if you inferred d=e=0 under DEC+J. Input a list of '' or NA to avoid this error.\n\n", 
                     sep = "")
    cat(errortxt)
    errortxt2 = paste("events_txt_list_into_events_table() is returning NULL which will might cause issues later.\n\n", 
                      sep = "")
    cat(errortxt2)
    return(NULL)
  }
  events_txt_list[is.na(events_txt_list)] = "none"
  noneTF = events_txt_list == "none"
  keepTF = (noneTF == FALSE)
  events_txt_list = events_txt_list[keepTF]
  if (length(events_txt_list) == 0) {
    events_table = NULL
    return(events_table)
  }
  if (length(trtable) > 0) {
    trtable_subset = NULL
  }
  tmptable = NULL
  for (i in 1:length(events_txt_list)) {
    tmptable_rows = events_txt_into_events_table(events_txt_list[i])
    rownums_in_trtable = as.numeric(tmptable_rows$nodenum_at_top_of_branch)
    num_newrows = nrow(tmptable_rows)
    tmptable = rbind(tmptable, tmptable_rows)
  }
  events_table = dfnums_to_numeric(adf2(tmptable))
  names(events_table) = c("nodenum_at_top_of_branch", 
                          "trynum", "brlen", "current_rangenum_1based", 
                          "new_rangenum_1based", "current_rangetxt", 
                          "new_rangetxt", "abs_event_time", "event_time", 
                          "event_type", "event_txt", "new_area_num_1based", 
                          "lost_area_num_1based", "dispersal_to", "extirpation_from")
  return(events_table)
}

events_txt_into_events_table<-function (branch_events_txt) 
{
  words = strsplit(branch_events_txt, split = ";")[[1]]
  events_table_for_branch = t(sapply(X = words, FUN = event_txt_into_events_row))
  row.names(events_table_for_branch) = NULL
  events_table_for_branch
  events_table_for_branch = adf2(events_table_for_branch)
  events_table_for_branch
  names(events_table_for_branch) = c("nodenum_at_top_of_branch", 
                                     "trynum", "brlen", "current_rangenum_1based", 
                                     "new_rangenum_1based", "current_rangetxt", 
                                     "new_rangetxt", "abs_event_time", "event_time", 
                                     "event_type", "event_txt", "new_area_num_1based", 
                                     "lost_area_num_1based", "dispersal_to", "extirpation_from")
  return(events_table_for_branch)
}

dfnums_to_numeric<-function (dtf, max_NAs = 0.5, printout = FALSE, roundval = NULL) 
{
  dtf_classes = cls.df(dtf, printout = FALSE)
  dtf_names = names(dtf)
  numcols = ncol(dtf)
  cls_col_list = c()
  for (i in 1:numcols) {
    cls_col = NA
    cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], 
                   "')", sep = "")
    eval(parse(text = cmdstr))
    cls_col_list[i] = cls_col
  }
  for (i in 1:numcols) {
    if (cls_col_list[i] == "list") {
      (next)()
    }
    if (cls_col_list[i] != "numeric") {
      newcol = NA
      cmdstr = paste("newcol = as.numeric(as.character(dtf$'", 
                     dtf_names[i], "'))", sep = "")
      suppressWarnings(eval(parse(text = cmdstr)))
      if (sum(is.na(newcol)) < (max_NAs * length(newcol))) {
        cmdstr = paste("dtf$'", dtf_names[i], "' = newcol", 
                       sep = "")
        suppressWarnings(eval(parse(text = cmdstr)))
        if (!is.null(roundval)) {
          cmdstr = paste("dtf$'", dtf_names[i], 
                         "' = round(dtf$'", dtf_names[i], "', digits=roundval)", 
                         sep = "")
          suppressWarnings(eval(parse(text = cmdstr)))
        }
      }
    }
  }
  tmp_classes = cls.df(dtf)
  dtf_classes$newclasses = tmp_classes[, ncol(tmp_classes)]
  if (printout) {
    cat("\n")
    cat("dfnums_to_numeric(dtf, max_NAs=", max_NAs, 
        ") reports: dataframe 'dtf_classes' has ", 
        nrow(dtf_classes), " rows, ", ncol(dtf_classes), 
        " columns.\n", sep = "")
    cat("...names() and classes() of each column below...\n", 
        sep = "")
    cat("\n")
    print(dtf_classes)
  }
  return(dtf)
}

adf2 <- function (x) 
{
  rownames = 1:nrow(x)
  return(as.data.frame(x, row.names = rownames, stringsAsFactors = FALSE))
}

event_txt_into_events_row<-function (word) 
{
  split_key_item <- function(word2) {
    output_pair = c("", "")
    words3 = strsplit(word2, split = ":")[[1]]
    numwords = length(words3)
    output_pair[1:numwords] = words3[1:numwords]
    return(output_pair)
  }
  words2 = strsplit(word, split = ",")[[1]]
  output = sapply(X = words2, FUN = split_key_item)
  tmprow = matrix(data = output[2, ], nrow = 1)
  return(tmprow)
}

cls.df <- function (dtf, printout = FALSE) 
{
  if (class(dtf) == "matrix") {
    dtf = as.data.frame(dtf, stringsAsFactors = FALSE)
  }
  dtf_names = names(dtf)
  numcols = ncol(dtf)
  cls_col_list = c()
  for (i in 1:numcols) {
    cls_col = NA
    cmdstr = paste("cls_col = class(dtf$'", dtf_names[i], 
                   "')", sep = "")
    eval(parse(text = cmdstr))
    cls_col_list[i] = cls_col
  }
  colnum = 1:numcols
  dtf_classes = cbind(colnum, dtf_names, cls_col_list)
  dtf_classes = data.frame(dtf_classes, row.names = colnum)
  if (printout) {
    cat("\n")
    cat("cls.df(dtf) reports: dataframe 'dtf' has ", 
        nrow(dtf), " rows, ", numcols, " columns.\n", 
        sep = "")
    cat("...names() and classes() of each column below...\n", 
        sep = "")
    cat("\n")
    print(dtf_classes)
    cat("\n")
  }
  return(dtf_classes)
}