#ifndef Chemistry_hpp
#define Chemistry_hpp

#include <string>
#include <vector>

namespace FL { namespace chem {

#ifdef __INTEL_COMPILER
    extern "C" {
        void chemistry_mp_avogadro_vibration_(const int & NAtoms, const char * symbol, const double * structure,
            const int & vibdim, const double * freq, const double * modeT, const char * file, const int & len_file);

        void chemistry_mp_initializephasefixing_(const int & NStates);
    }

    inline void Avogadro_Vibration(const int & NAtoms, const std::vector<std::string> & symbol, const double * structure,
        const int & vibdim, const double * freq, const double * modeT, const std::string & file="avogadro.log"
    ) {
        char * s = new char[2 * symbol.size()];
        for (size_t i = 0; i < symbol.size(); i++) {
            s[2 * i] = symbol[i].c_str()[0];
            if (symbol[i].size() > 1) s[2 * i + 1] = symbol[i].c_str()[1];
            else                      s[2 * i + 1] = ' ';
        }
        chemistry_mp_avogadro_vibration_(NAtoms, s, structure, vibdim, freq, modeT, file.c_str(), file.size());
        delete [] s;
    }

    inline void InitializePhaseFixing(const int & NStates) {chemistry_mp_initializephasefixing_(NStates);}
#elif __GNUC__
    extern "C" {
        void __chemistry_MOD_avogadro_vibration(const int & NAtoms, const char * symbol, const double * structure,
            const int & vibdim, const double * freq, const double * modeT, const char * file, const int & len_file);

        void __chemistry_MOD_initializephasefixing(const int & NStates);
    }

    inline void Avogadro_Vibration(const int & NAtoms, const std::vector<std::string> & symbol, const double * structure,
        const int & vibdim, const double * freq, const double * modeT, const std::string & file="avogadro.log"
    ) {
        char * s = new char[2 * symbol.size()];
        for (size_t i = 0; i < symbol.size(); i++) {
            s[2 * i] = symbol[i].c_str()[0];
            if (symbol[i].size() > 1) s[2 * i + 1] = symbol[i].c_str()[1];
            else                      s[2 * i + 1] = ' ';
        }
        __chemistry_MOD_avogadro_vibration(NAtoms, s, structure, vibdim, freq, modeT, file.c_str(), file.size());
        delete [] s;
    }

    inline void InitializePhaseFixing(const int & NStates) {__chemistry_MOD_initializephasefixing(NStates);}
#endif

} // namespace chem
} // namespace FL

#endif