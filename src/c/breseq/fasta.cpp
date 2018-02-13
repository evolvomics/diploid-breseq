/*****************************************************************************
 
 AUTHORS
 
 Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
 David B. Knoester
 
 LICENSE AND COPYRIGHT
 
 Copyright (c) 2008-2010 Michigan State University
 Copyright (c) 2011-2017 The University of Texas at Austin
 
 breseq is free software; you can redistribute it and/or modify it under the
 terms the GNU General Public License as published by the Free Software
 Foundation; either version 1, or (at your option) any later version.
 
 *****************************************************************************/

#include "libbreseq/fasta.h"

using namespace std;

namespace breseq {
  
  //constructor
  cFastaFile::cFastaFile(const string &file_name, ios_base::openmode mode) 
    : fstream(file_name.c_str(), mode)
    , m_file_name(file_name)
    , m_current_line()
    , m_current_line_num(0)
    , m_bases_per_line(60)
  {
    ASSERT(!(*this).fail(), "Failed to open FASTA file: " + m_file_name);
    
    if (mode == ios_base::in) {
      breseq::getline(*this, m_current_line);
    }
  }

  // read one sequence record from the file
  bool cFastaFile::read_sequence(cFastaSequence &sequence) {
    
    // clear sequence
    sequence.set_name("");
    sequence.set_description("");
    sequence.set_sequence("");
    
    // We're done, no error
    if (this->eof()) return false;
    
    // Current line should begin with >
    assert(m_current_line[0] == '>');
    
    // @JEB it may be better to truncate at the space
    
    // The sequence name is the first word
    size_t pos = m_current_line.find_first_of(" \t\r\n", 1);
    sequence.set_name(m_current_line.substr(1,(pos != string::npos) ? pos-1 : string::npos));
    pos = m_current_line.find_first_not_of( " \t\r\n", pos);
    if (pos != string::npos) sequence.set_description(m_current_line.substr(pos));
        
    breseq::getline(*this, m_current_line);
    m_current_line_num++;
    string nucleotide_sequence;
    while ((m_current_line[0] != '>') && !this->eof()) {
      
      // Clean the sequence of spaces and extra line returns ('\r' is particularly dangerous).
      // We could also check for valid characters
      
      m_current_line = substitute(m_current_line, "\r", "");
      m_current_line = substitute(m_current_line, "\n", "");
      
      
      nucleotide_sequence += m_current_line;

      breseq::getline(*this, m_current_line);
      m_current_line_num++;
    }
    
    assert(nucleotide_sequence.length() > 0);
    sequence.set_sequence(nucleotide_sequence);

    return true;
  }

  void cFastaFile::write_sequence(const cFastaSequence &sequence, const int32_t chr_0) {
    
    // If you request a certain chromosome, we suffix the name with that number
    int32_t chr_to_print = chr_0;
    if (chr_to_print == -1) {
      (*this) << ">" << sequence.get_name() << endl;
      chr_to_print=0;
    } else {
      ASSERT( (chr_0 >= 0) && (static_cast<size_t>(chr_0) < sequence.get_ploidy()), "Requested chromosome number (" + to_string(chr_0+1) + ") that does not exist. Ploidy of sequence is " + to_string(sequence.get_ploidy()));
      (*this) << ">" << sequence.get_name() << ":" << chr_0 << endl;
    }
    
    int32_t start_1 = 1;
    while (start_1 <= (int32_t)sequence.get_sequence_length(chr_to_print)) {
      (*this) << sequence.get_sequence_1(start_1, min<int32_t>(start_1 + m_bases_per_line - 1, sequence.get_sequence_length()), chr_to_print) << endl;
      start_1 += m_bases_per_line;
    }

  }
  
  
} // breseq namespace

