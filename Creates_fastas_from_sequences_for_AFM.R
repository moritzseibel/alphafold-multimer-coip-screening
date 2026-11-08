#get domain sequence from "full_seq"
extract_domain <- function(full_seq, range) {

  substr(full_seq, range[1], range[2])
}

#writes fasta file from sequences A (NPL4) and B (OTUD5)
write_fasta_AB <- function(domain_seq, peptide_seq, file) {
  con <- file(file, "w")
  writeLines(">A", con)
  writeLines(domain_seq, con)
  writeLines(">B", con)
  writeLines(peptide_seq, con)
  close(con)
}

#define the protein sequences manually (fetch from pdb would be possible, but unnecessary)
OTUD5Seq <- "MTILPKKKPPPPDADPANEPPPPGPMPPAPRRGGGVGVGGGGTGVGGGDRDRDSGVVGARPRASPPPQGPLPGPPGALHRWALAVPPGAVAGPRPQQASPPPCGGPGGPGGGPGDALGAAAAGVGAAGVVVGVGGAVGVGGCCSGPGHSKRRRQAPGVGAVGGGSPEREEVGAGYNSEDEYEAAAARIEAMDPATVEQQEHWFEKALRDKKGFIIKQMKEDGACLFRAVADQVYGDQDMHEVVRKHCMDYLMKNADYFSNYVTEDFTTYINRKRKNNCHGNHIEMQAMAEMYNRPVEVYQYSTGTSAVEPINTFHGIHQNEDEPIRVSYHRNIHYNSVVNPNKATIGVGLGLPSFKPGFAEQSLMKNAIKTSEESWIEQQMLEDKKRATDWEATNEAIEEQVARESYLQWLRDQEKQARQVRGPSQPRKASATCSSATAAASSGLEEWTSRSPRQRSSASSPEHPELHAELGMKPPSPGTVLALAKPPSPCAPGTSSQFSAGADRATSPLVSLYPALECRALIQQMSPSAFGLNDWDDDEILASVLAVSQQEYLDSMKKNKVHRDPPPDKS"

NPLOC4Seq <- "MAESIIIRVQSPDGVKRITATKRETAATFLKKVAKEFGFQNNGFSVYINRNKTGEITASSNKSLNLLKIKHGDLLFLFPSSLAGPSSEMETSVPPGFKVFGAPNVVEDEIDQYLSKQDGKIYRSRDPQLCRHGPLGKCVHCVPLEPFDEDYLNHLEPPVKHMSFHAYIRKLTGGADKGKFVALENISCKIKSGCEGHLPWPNGICTKCQPSAITLNRQKYRHVDNIMFENHTVADRFLDFWRKTGNQHFGYLYGRYTEHKDIPLGIRAEVAAIYEPPQIGTQNSLELLEDPKAEVVDEIAAKLGLRKVGWIFTDLVSEDTRKGTVRYSRNKDTYFLSSEECITAGDFQNKHPNMCRLSPDGHFGSKFVTAVATGGPDNQVHFEGYQVSNQCMALVRDECLLPCKDAPELGYAKESSSEQYVPDVFYKDVDKFGNEITQLARPLPVEYLIIDITTTFPKDPVYTFSISQNPFPIENRDVLGETQDFHSLATYLSQNTSSVFLDTISDFHLLLFLVTNEVMPLQDSISLLLEAVRTRNEELAQTWKRSEQWATIEQLCSTVGGQLPGLHEYGAVGGSTHTATAAMWACQHCTFMNQPGTGHCEMCSLPRT"

#define the domain sections manually
OTUD5Domains <- list(domain1 = c(173, 346), domain2 = c(364, 419), domain3 = c(539,559))
NPLOC4Domains <- list(domain1 = c(3,80), domain2 = c(104, 559), domain3= c(582, 608))

#cuts away the name of the domain names, to leave only the coordinates
OTUD5DomVec <- unname(OTUD5Domains)
NPLOC4DomVec <- unname(NPLOC4Domains)

#split peptides in len = 20 and step = 5
split_peptides <- function(seq, len = 20, step = 5, domains = NULL) {
  peptides <- list()
  num <- 1
  i <- 1
  seq_len <- nchar(seq)
  
  while (i + len - 1 <= seq_len) {
    skip <- FALSE
    
#if start would be inside domain -> jump to 1. aa after domain
    if (!is.null(domains)) {
      for (d in domains) {
        if (i >= d[1] && i <= d[2]) {
          i <- d[2] + 1
          skip <- TRUE
          break
        }
      }
    }
    if (skip) next
    
    end <- i + len - 1
#defines, where we should jump after peptide
    jump_after <- NULL
    
#check, if peptide would "contain" a domain part
    if (!is.null(domains)) {
      for (d in domains) {
        if (end >= d[1] && i < d[1]) {
#cut peptide just before domain
          end <- d[1] - 1
#next peptide starts after
          jump_after <- d[2] + 1
          break
        }
      }
    }
    
    peptides[[paste0("Peptide_", num)]] <- list(
      seq   = substr(seq, i, end),
      start = i
    )
    num <- num + 1
    
#continue after domain, if recent peptide was cut due to domain presence
    if (!is.null(jump_after)) {
      i <- jump_after
    } else {
      i <- i + step
    }
  }
  
  peptides
}


peptides_OTUD5 <- split_peptides(OTUD5Seq, len = 20, step = 5, domains = OTUD5DomVec)
peptides_NPLOC4 <- split_peptides(NPLOC4Seq, len = 20, step = 5, domains = NPLOC4DomVec)


#creates folder for fasta output
output_dir <- "Scripts_output/fasta_for_2nd_run"
if (!dir.exists(output_dir)) dir.create(output_dir)

#Loops all domain fasta pairs
#NPLOC4-Domänen + OTUD5-Peptide
for (d_idx in seq_along(NPLOC4DomVec)) {
  domain_range <- NPLOC4DomVec[[d_idx]]
  domain_seq <- extract_domain(NPLOC4Seq, domain_range)
  domain_start <- domain_range[1]
  for (p_idx in seq_along(peptides_OTUD5)) {
    peptide_info <- peptides_OTUD5[[p_idx]]
    peptide_seq <- peptide_info$seq
    peptide_start <- peptide_info$start
#domain boundaries
    dom_start <- domain_range[1]
    dom_end <- domain_range[2]
    dom_start_aa <- substr(NPLOC4Seq, dom_start, dom_start)
    dom_end_aa <- substr(NPLOC4Seq, dom_end, dom_end)
#peptide boundaries
    pep_start <- peptide_info$start
    pep_end <- peptide_info$start + nchar(peptide_seq) - 1
    pep_start_aa <- substr(OTUD5Seq, pep_start, pep_start)
    pep_end_aa <- substr(OTUD5Seq, pep_end, pep_end)
#Filename
    job_name <- paste0("NPLOC4_domain_", d_idx,"(",dom_start_aa,dom_start,"-",dom_end_aa,dom_end, ")_OTUD5_peptide", p_idx,"(",pep_start_aa,pep_start,"-",pep_end_aa,pep_end,")")
    file_path <- file.path(output_dir, paste0(job_name, ".fasta"))
#creates fastas
    write_fasta_AB(
      domain_seq  = domain_seq,
      peptide_seq = peptide_seq,
      file        = file_path
    )
    
  }
}

#OTUD5-Domänen + NPLOC4-Peptide
for (d_idx in seq_along(OTUD5DomVec)) {
  domain_range <- OTUD5DomVec[[d_idx]]
  domain_seq <- extract_domain(OTUD5Seq, domain_range)
  domain_start <- domain_range[1]
  for (p_idx in seq_along(peptides_NPLOC4)) {
    peptide_info <- peptides_NPLOC4[[p_idx]]
    peptide_seq <- peptide_info$seq
    peptide_start <- peptide_info$start
    
#Domain boundaries
    dom_start <- domain_range[1]
    dom_end   <- domain_range[2]
    dom_start_aa <- substr(OTUD5Seq, dom_start, dom_start)
    dom_end_aa   <- substr(OTUD5Seq, dom_end,   dom_end)
  
#Peptide boundaries
    pep_start <- peptide_info$start
    pep_end   <- peptide_info$start + nchar(peptide_seq) - 1
    pep_start_aa <- substr(NPLOC4Seq, pep_start, pep_start)
    pep_end_aa   <- substr(NPLOC4Seq, pep_end,   pep_end)
    
#Filename
    job_name <- paste0("OTUD5_domain_", d_idx,"(",dom_start_aa,dom_start,"-",dom_end_aa,dom_end, ")_NPLOC4_peptide", p_idx,"(",pep_start_aa,pep_start,"-",pep_end_aa,pep_end,")")
    file_path <- file.path(output_dir, paste0(job_name, ".fasta"))
    
#create fasta
    write_fasta_AB(
      domain_seq  = domain_seq,
      peptide_seq = peptide_seq,
      file        = file_path
    )
    
  }
}

#Full-length NPLOC4 vs full-length OTUD5 
job_name  <- "NPLOC4_full_M1-T608_vs_OTUD5_full_M1-S571"
file_path <- file.path(output_dir, paste0(job_name, ".fasta"))

#create fasta
write_fasta_AB(
  domain_seq  = NPLOC4Seq,  # chain A
  peptide_seq = OTUD5Seq,   # chain B
  file        = file_path
)

#OTUD5 domains x NPLOC4 domains
for (d_idx in seq_along(OTUD5DomVec)) {
  domain_range <- OTUD5DomVec[[d_idx]]
  domain_seq <- extract_domain(OTUD5Seq, domain_range)
  domain_start <- domain_range[1]
  for (n_idx in seq_along(NPLOC4DomVec)) {
    domain_range2 <- NPLOC4DomVec[[n_idx]]
    domain_seq2 <- extract_domain(NPLOC4Seq, domain_range2)
    domain_start2 <- domain_range2[1]
    
#Domain boundaries
    dom_start <- domain_range[1]
    dom_end   <- domain_range[2]
    dom_start_aa <- substr(OTUD5Seq, dom_start, dom_start)
    dom_end_aa   <- substr(OTUD5Seq, dom_end,   dom_end)
    
#Domain boundaries
    dom_start2 <- domain_range2[1]
    dom_end2   <- domain_range2[2]
    dom_start_aa2 <- substr(NPLOC4Seq, dom_start2, dom_start2)
    dom_end_aa2   <- substr(NPLOC4Seq, dom_end2,   dom_end2)
    
#Filename
    job_name <- paste0("OTUD5_domain_", d_idx,"(",dom_start_aa,dom_start,"-",dom_end_aa,dom_end, ")_NPLOC4_domain", n_idx,"(",dom_start_aa2,dom_start2,"-",dom_end_aa2,dom_end2,")")
    file_path <- file.path(output_dir, paste0(job_name, ".fasta"))
    
#create fasta
    write_fasta_AB(
      domain_seq  = domain_seq,
      peptide_seq = domain_seq2,
      file        = file_path
    )
    
  }
}

#OTUD5 peptides & domains x NPLOC4 whole
np_full_start    <- 1
np_full_end      <- nchar(NPLOC4Seq)
np_full_start_aa <- substr(NPLOC4Seq, np_full_start, np_full_start)
np_full_end_aa   <- substr(NPLOC4Seq, np_full_end,   np_full_end)

#NPLOC4 full x OTUD5 peptides 
for (p_idx in seq_along(peptides_OTUD5)) {
  peptide_info <- peptides_OTUD5[[p_idx]]
  peptide_seq  <- peptide_info$seq
  
  pep_start    <- peptide_info$start
  pep_end      <- peptide_info$start + nchar(peptide_seq) - 1
  pep_start_aa <- substr(OTUD5Seq, pep_start, pep_start)
  pep_end_aa   <- substr(OTUD5Seq, pep_end,   pep_end)
  
  job_name <- paste0(
    "NPLOC4_full(",
    np_full_start_aa, np_full_start, "-",
    np_full_end_aa,   np_full_end, ")",
    "_OTUD5_peptide", p_idx, "(",
    pep_start_aa, pep_start, "-",
    pep_end_aa,   pep_end, ")"
  )
  
  file_path <- file.path(output_dir, paste0(job_name, ".fasta"))
  write_fasta_AB(
    domain_seq  = NPLOC4Seq,   #chain A = full NPLOC4
    peptide_seq = peptide_seq, #chain B = OTUD5 peptide
    file        = file_path
  )
}

#NPLOC4 full x OTUD5 domains 
for (d_idx in seq_along(OTUD5DomVec)) {
  domain_range <- OTUD5DomVec[[d_idx]]
  domain_seq   <- extract_domain(OTUD5Seq, domain_range)
  
  dom_start    <- domain_range[1]
  dom_end      <- domain_range[2]
  dom_start_aa <- substr(OTUD5Seq, dom_start, dom_start)
  dom_end_aa   <- substr(OTUD5Seq, dom_end,   dom_end)
  
  job_name <- paste0(
    "NPLOC4_full(",
    np_full_start_aa, np_full_start, "-",
    np_full_end_aa,   np_full_end, ")",
    "_OTUD5_domain_", d_idx, "(",
    dom_start_aa, dom_start, "-",
    dom_end_aa,   dom_end, ")"
  )
  
  file_path <- file.path(output_dir, paste0(job_name, ".fasta"))
  write_fasta_AB(
    domain_seq  = NPLOC4Seq,  #chain A = full NPLOC4
    peptide_seq = domain_seq, #chain B = OTUD5 domain
    file        = file_path
  )
}