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

#ifndef _BRESEQ_FASTA_H_
#define _BRESEQ_FASTA_H_

#include "common.h"

using namespace std;

namespace breseq {
	
	/*! Interface for loading and manipulating files in FASTA format.
   */
  
  /*! Sequence class.
   */
   
  class cFastaSequence {
    
    private:
      string m_name;                //>NAME DESCRIPTION
      string m_description;         //
      vector<string> m_sequence;    //sequence ... multiple because we use this to store multiploid reference genome info
    
    public:
      cFastaSequence() {}
      cFastaSequence(const string& _name, const string& _description, const string& _sequence)
    : m_name(_name), m_description(_description) { m_sequence.clear(); m_sequence.push_back(_sequence); }
      ~cFastaSequence() {}
    
    void set_name(const string& _name) { m_name = _name; }
    void set_description(const string& _description) { m_description = _description; }
    void set_sequence(const string& _sequence, const size_t chr_0=0)
    {
      if (chr_0 + 1 > m_sequence.size()) m_sequence.resize(chr_0+1,"");
      m_sequence[chr_0] = _sequence;
    }
    
    const string get_description() const { return m_description; }
    const string get_name() const { return m_name; }
    const string get_sequence(const size_t _chromosome_index=0) const { return m_sequence[_chromosome_index]; }
    
    // Utility function for checking and correcting bounds
    void correct_coordinate_bounds(int64_t &pos_1, const size_t chr_0) const {
      if (pos_1 < 1) {
        WARN_WITH_BACKTRACE("Coordinate (" + to_string<int64_t>(pos_1) + ") requested for get_sequence_1 is < 1.");
        pos_1 = 1;
      }
      
      if (pos_1 > (int64_t)(m_sequence[chr_0].length())) {
        WARN_WITH_BACKTRACE("Coordinate (" + to_string<int64_t>(pos_1) + ") requested for get_sequence_1 is > length of sequence (" + to_string(m_sequence[chr_0].length()) + ") .");
        pos_1 = m_sequence[chr_0].length();
      }
      
    }
    
    void correct_range_bounds(int64_t &start_1, int64_t &end_1) const {
      
      if (start_1 > end_1) {
        WARN_WITH_BACKTRACE("Start coordinate (" + to_string<int64_t>(start_1) + ") is > end coordinate (" + to_string<int64_t>(end_1) + ") for range.");
        end_1 = start_1;
      }
    }
    
    // We do some checking here to make sure we don't throw an out-of-bounds error.
    const string get_sequence_1(int64_t start_1, int64_t end_1, const size_t chr_0=0) const
    {
      
      // No error, this is sometimes intended
      if ((start_1==0) && (end_1==0)) return "";
      
      correct_coordinate_bounds(start_1, chr_0);
      correct_coordinate_bounds(end_1, chr_0);
      correct_range_bounds(start_1, end_1);
      return m_sequence[chr_0].substr(start_1-1, end_1-start_1+1);
    }
    
    char get_char_1(int64_t pos_1, const size_t chr_0=0) const
    {
      correct_coordinate_bounds(pos_1, chr_0);
      return m_sequence[chr_0][pos_1-1];
    }

    
    void replace_sequence_1(int64_t start_1, int64_t end_1, const string &replacement_seq, const size_t chr_0=0)
    {
      correct_coordinate_bounds(start_1, chr_0);
      correct_coordinate_bounds(end_1, chr_0);
      correct_range_bounds(start_1, end_1);
      m_sequence[chr_0].replace(start_1-1, end_1-start_1+1, replacement_seq);
    }
    
    void insert_sequence_1(int64_t pos_1, const string &insertion_seq, const size_t chr_0=0)
    {
      // Allow a value of zero, which means insert at the start of the sequence
      if (pos_1 != 0) correct_coordinate_bounds(pos_1, chr_0);
      m_sequence[chr_0].insert(pos_1, insertion_seq);
    }

    
    size_t get_sequence_length(const size_t chr_0=0) const
    {
      return m_sequence[chr_0].length();
    }
    
    void set_ploidy(const size_t _ploidy)
    {
      for(size_t i=1;i<_ploidy;i++) {
        m_sequence.push_back(m_sequence[0]);
      }
    }
    
    size_t get_ploidy() const
    {
      return m_sequence.size();
    }
    
    void standardize_sequence();

   };
  
  inline ostream& operator<<(ostream& out, const cFastaSequence& fasta_sequence)
  {
    out << ">" <<fasta_sequence.get_name() << endl;
    out << fasta_sequence.get_sequence() << endl;
    return out;
  }


	/*! File class.
	 */ 
  
  class cFastaFile : public fstream {
    
  public:
    string    m_file_name;  
    
  protected:
    string    m_current_line;
    uint32_t  m_current_line_num;
    uint32_t  m_bases_per_line;
    
  public:  
    cFastaFile(const string &file_name, ios_base::openmode mode);
    ~cFastaFile() {};
      
    bool read_sequence(cFastaSequence &sequence);
    void write_sequence(const cFastaSequence &sequence, const int32_t chr_0=-1);
    void set_current_line(const string& current_line)
    { m_current_line = current_line;}

  
  };
	
} // breseq namespace

#endif
