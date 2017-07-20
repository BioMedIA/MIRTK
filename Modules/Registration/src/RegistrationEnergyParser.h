/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mirtk/String.h"
#include "mirtk/Stream.h"
#include "mirtk/Array.h"


namespace mirtk {


/**
 * GRAMMAR (not up to date and complete)
 *
 * energy:
 *   similarities
 *   similarities + constraints
 *
 * similarities:
 *  weighted_similarity + similarities
 *  weighted_similarity
 *
 * weighted_similiarity:
 *   weight * similarity
 *   weight similarity
 *   similarity
 *
 * similarity:
 *   NAME(NAME, NAME)
 *   NAME[NAME](NAME, NAME)
 *
 * constraints:
 *   weighted_constraint + constraints
 *   weighted_constraint
 *
 * weighted_constraint:
 *   weight * constraint
 *   weight constraint
 *   constraint
 *
 * constraint:
 *   NAME
 *   NAME(T)
 *   NAME[NAME]
 *   NAME[NAME](T)
 *
 * weight:
 *   NUMBER
 *   NUMBER/NUMBER
 *
 */
class RegistrationEnergyParser
{
  GenericRegistrationFilter *_Filter;

  typedef GenericRegistrationFilter::TransformationInfo     TransformationInfo;
  typedef GenericRegistrationFilter::ImageSimilarityInfo    ImageSimilarityInfo;
  typedef GenericRegistrationFilter::PointSetDistanceInfo   PointSetDistanceInfo;
  typedef GenericRegistrationFilter::PointSetConstraintInfo PointSetConstraintInfo;
  typedef GenericRegistrationFilter::ConstraintInfo         ConstraintInfo;

  // ---------------------------------------------------------------------------
  /// Enumeration of input data types
  enum DataType
  {
    IMAGE,
    POLYDATA
  };

  // ---------------------------------------------------------------------------
  /// Enumeration of terminal symbols
  enum TokenEnum
  {
    NAME,
    INTEGER,
    NUMBER,
    END,
    SPACE  = ' ',
    PLUS   = '+',
    MINUS  = '-',
    MUL    = '*',
    DIV    = '/',
    LP     = '(',
    RP     = ')',
    LB     = '[',
    RB     = ']',
    LCB    = '{',
    RCB    = '}',
    COMMA  = ',',
    COLON  = ':',
    CIRC   = 'o', ///< Function composition
    CARET  = '^'  ///< Exponentiation of transformation
  };

  // ---------------------------------------------------------------------------
  /// Token of energy function description
  struct Token
  {
    TokenEnum   _Token;
    double      _Number;
    int         _Integer;
    string _Name;

    Token(TokenEnum token)  : _Token(token)                      {}
    Token(char      c)      : _Token(TokenEnum(c))               {}
    Token(double    number) : _Token(NUMBER), _Number(number),
                              _Integer(static_cast<int>(number)) {}
    Token(string name) : _Token(NAME),   _Name(name)        {}

    Token(const Token &token)
    :
      _Token  (token._Token),
      _Number (token._Number),
      _Integer(token._Integer),
      _Name   (token._Name)
    {}

    Token &operator =(const Token &token)
    {
      _Token   = token._Token;
      _Number  = token._Number;
      _Integer = token._Integer;
      _Name    = token._Name;
      return *this;
    }

    bool operator ==(TokenEnum token) const { return _Token == token; }
    bool operator !=(TokenEnum token) const { return _Token != token; }
  };

  // ---------------------------------------------------------------------------
  /// Get next token of energy function description
  Token NextToken(istream &in, bool skip_ws = true)
  {
    char c;

    // Skip whitespaces and detect end of input
    if (skip_ws) {
      do {
        if (!in.get(c)) return END;
      } while (::isspace(c));
    } else {
      if (!in.get(c)) return END;
    }

    // Get next token
    switch (c) {
      case ' ':
      case '*':
      case '/':
      case '-':
      case '+':
      case '(':
      case ')':
      case '[':
      case ']':
      case ',':
      case ':':
      case '^':
      case '{':
      case '}':
        return c;
      case '0': case '1': case '2': case '3': case '4':
      case '5': case '6': case '7': case '8': case '9':
      case '.': {
        string number(1, c);
        bool before_decimal_point = (c != '.');
        while (in.get(c)) {
          if (isdigit(c) || (c == '.' && before_decimal_point)) {
            if (c == '.') before_decimal_point = false;
            number += c;
          } else if (c == 'e' || c == 'E') {
            bool is_next_token = false;
            string exponent(1, c);
            if (in.get(c)) {
              exponent += c;
              if (c == '-' || c == '+') {
                if (in.get(c)) {
                  exponent += c;
                  if (!isdigit(c)) {
                    is_next_token = true;
                  }
                } else {
                  is_next_token = true;
                }
              } else if (!isdigit(c)) {
                is_next_token = true;
              }
              if (!is_next_token) {
                while (in.get(c)) {
                  if (isdigit(c)) {
                    exponent += c;
                  } else {
                    in.putback(c);
                    break;
                  }
                }
              }
            } else {
              is_next_token = true;
            }
            if (is_next_token) {
              for (auto it = exponent.rbegin(); it != exponent.rend(); ++it) {
                in.putback(*it);
              }
            } else {
              number += exponent;
            }
            break;
          } else {
            in.putback(c);
            break;
          }
        }
        double value = .0;
        if (!FromString(number, value)) {
          in.setstate(std::ios::failbit);
        }
        return value;
      }
      case 'o': {
        char c2;
        if (!in.get(c2)) return c;
        in.putback(c2);
        if (!::isalnum(c2) && c2 != '_') return c;
        string name(1, c);
        while (in.get(c) && (isalnum(c) || c == '_')) name += c;
        in.putback(c);
        return name;
      }
      default: {
        if (::isalpha(c)) {
          string name(1, c);
          while (in.get(c) && (isalnum(c) || c == '_')) name += c;
          in.putback(c);
          return name;
        }
        cout << endl;
        cerr << "Invalid character in objective function description: \"" << c
                  << "\" (ASCII code " << int(c) << ")" << endl;
        exit(1);
        return END;
      }
    }
  }

  // ---------------------------------------------------------------------------
  /// Evaluate weight of energy term
  double Number(istream &in, Token &token)
  {
    double w = 1.0;
    while (token == MINUS || token == PLUS) {
      if (token == MINUS) w *= -1.0;
      token = NextToken(in);
    }
    if (token == NUMBER) {
      w *= token._Number;
      token = NextToken(in);
      if (token != DIV) return w;
      token = NextToken(in);
      if (token != NUMBER) {
        cout << endl;
        cerr << "Expected number after / in energy formula" << endl;
        exit(1);
      }
      double d = token._Number;
      if (d == .0) {
        cout << endl;
        cerr << "Division by zero in energy formula!" << endl;
        exit(1);
      }
      token = NextToken(in);
      return w / d;
    }
    return w;
  }

  // ---------------------------------------------------------------------------
  /// Whether a given number is an integer value
  bool IsInteger(double number)
  {
    return (static_cast<double>(static_cast<int>(number)) == number);
  }

  // ---------------------------------------------------------------------------
  /// Parse input identifier with MATLAB style indexing syntax
  Array<int> InputIndex(istream &in, Token &token, const string &name,
                        DataType type, int max_index)
  {
    Array<int> index;
    const char *c;
    const char *type_name;
    const char *id_name;

    switch (type) {
      case IMAGE:
        c         = "I"; // as in "Image"
        type_name = "images";
        id_name   = "image identifier";
        break;
      case POLYDATA:
        c         = "PCS"; // as in "Point set", "Curve", "Surface"
        type_name = "point set";
        id_name   = "point set identifier";
        break;
      default:
        cerr << "RegistrationEnergyParser::InputIndex: Invalid data type identifier: " << type << endl;
        exit(1);
    }

    if (token != NAME || token._Name.substr(0, 1).find_first_of(c) != 0) return index;

    if (token._Name.size() == 1) {
      token = NextToken(in);
      if (token != LP) {
        cout << endl;
        cerr << "Expected single index or ( after input identifier of similarity term: " << name << endl;
        exit(1);
      }
      token = NextToken(in);
      bool inside_square_brackets = (token == LB);
      if (inside_square_brackets) token = NextToken(in);
      do {
        int start = 1;
        if (token == NAME && token._Name == "end") {
          start = -1;
          token = NextToken(in);
          if (token == MINUS) {
            token = NextToken(in);
            if (token != NUMBER || !IsInteger(token._Number) || token._Integer < 0) {
              cout << endl;
              cerr << "Expected positive integer after 'end-' in input identifier of similarity term: " << name << endl;
              exit(1);
            }
            start = -token._Integer;
            token = NextToken(in);
          }
        } else if (token == NUMBER && IsInteger(token._Number)) {
          start = token._Integer;
          token = NextToken(in);
        } else {
          cout << endl;
          cerr << "Expected input index or 'end' in identifier argument list of similarity term: " << name << endl;
          exit(1);
        }
        int inc = 1;
        int end = start;
        if (token == COLON) {
          token = NextToken(in);
          if (token == NAME && token._Name == "end") {
            end = -1;
            token = NextToken(in);
            if (token == MINUS) {
              token = NextToken(in);
              if (token != NUMBER || !IsInteger(token._Number) || token._Integer < 0) {
                cout << endl;
                cerr << "Expected positive integer after 'end-' in input identifier of similarity term: " << name << endl;
                exit(1);
              }
              end   = -token._Integer;
              token = NextToken(in);
            }
          } else if (token == NUMBER) {
            if (!IsInteger(token._Number) || token._Integer < 0) {
              cout << endl;
              cerr << "Expected positive integer after : in input identifier of similarity term: " << name << endl;
              exit(1);
            }
            end = token._Integer;
            token = NextToken(in);
          } else {
            cout << endl;
            cerr << "Expected input index, 'end' or increment after : in input identifier of similarity term: " << name << endl;
            exit(1);
          }
          if (token == COLON) {
            inc   = end;
            end   = start;
            token = NextToken(in);
            if (token == NAME && token._Name == "end") {
              end = -1;
              token = NextToken(in);
              if (token == MINUS) {
                token = NextToken(in);
                if (token != NUMBER || !IsInteger(token._Number) || token._Integer < 0) {
                  cout << endl;
                  cerr << "Expected positive integer after 'end-' in input identifier of similarity term: " << name << endl;
                  exit(1);
                }
                end   = -token._Integer;
                token = NextToken(in);
              }
            } else if (token == NUMBER) {
              if (!IsInteger(token._Number) || token._Integer < 0) {
                cout << endl;
                cerr << "Expected positive integer after : in input identifier of similarity term: " << name << endl;
                exit(1);
              }
              end   = token._Integer;
              token = NextToken(in);
            } else {
              cout << endl;
              cerr << "Expected input index or 'end' after second : in input identifier of similarity term: " << name << endl;
              exit(1);
            }
          }
        }
        if (end < 0) end += max_index + 1;
        if (start < 1 || end < start || inc == 0) {
          cout << endl;
          cerr << "Invalid index (range) in input identifier of energy term: " << name << endl;
          exit(1);
        }
        if (start > max_index || end > max_index) {
          cout << endl;
          cerr << "Not enough input " << type_name << " or invalid index (range) in " << id_name << " of similarity term: " << name << endl;
          exit(1);
        }
        for (int i = start - 1; i < end; i += inc) index.push_back(i);
        if (token == RB) {
          if (!inside_square_brackets) {
            cout << endl;
            cerr << "Found ] without matching [ in input identifier of similarity term: " << name << endl;
            exit(1);
          }
          inside_square_brackets = false;
          token = NextToken(in);
        } else if (token == END) {
          cout << endl;
          cerr << "Expected ] after input index (range) in input identifier of similarity term: " << name << endl;
          exit(1);
        }
      } while (inside_square_brackets);
      if (token != RP) {
        cout << endl;
        cerr << "Expected ) after input index (range) in input identifier of similarity term: " << name << endl;
        exit(1);
      }
      token = NextToken(in);
    } else {
      int i = -1;
      if (FromString(token._Name.c_str() + 1, i) && i > 0 &&
          token._Name == (token._Name.substr(0, 1) + ToString(i))) {
        if (i > max_index) {
          cout << endl;
          cerr << "Not enought input " << type_name << " or invalid index in " << id_name << " of similarity term: " << name << endl;
          exit(1);
        }
        index.push_back(i - 1);
        token = NextToken(in);
      } else {
        cout << endl;
        cerr << "Not enough input " << type_name << " or invalid " << id_name << " " << token._Name << " in argument list of similarity term: " << name << endl;
        exit(1);
      }
    }
    return index;
  }

public:

  // ---------------------------------------------------------------------------
  /// Substitute single substring by given value
  template <class T>
  static string Substitute(const string &s, const char *var, T value)
  {
    size_t pos = s.find(var);
    if (pos == string::npos) return s;
    return s.substr(0, pos) + ToString(value) + s.substr(pos + strlen(var));
  }

protected:

  // ---------------------------------------------------------------------------
  /// Name of image dissimilarity term
  string TermName(const string &str, const ImageSimilarityInfo &info, int i = -1) const
  {
    string name(str);
    name = Substitute(name, "{t}", info._TargetIndex + 1);
    name = Substitute(name, "{s}", info._SourceIndex + 1);
    if (i > -1) name = Substitute(name, "{i}", i + 1);
    return name;
  }

  // ---------------------------------------------------------------------------
  /// Name of point set distance term
  string TermName(const string &str, const PointSetDistanceInfo &info, int i = -1) const
  {
    string name(str);
    name = Substitute(name, "{t}", info._TargetIndex + 1);
    name = Substitute(name, "{s}", info._SourceIndex + 1);
    if (i > -1) name = Substitute(name, "{i}", i + 1);
    return name;
  }

  // ---------------------------------------------------------------------------
  /// Name of point set constraint term
  string TermName(const string &str, const PointSetConstraintInfo &info, int i = -1) const
  {
    string name(str);
    name = Substitute(name, "{n}", info._PointSetIndex + 1);
    name = Substitute(name, "{t}", info._PointSetIndex + 1);
    if (info._RefPointSetIndex > -1) {
      name = Substitute(name, "{s}", info._RefPointSetIndex + 1);
    } else if (info._RefImageIndex > -1) {
      name = Substitute(name, "{s}", info._RefImageIndex + 1);
    }
    if (i > -1) name = Substitute(name, "{i}", i + 1);
    return name;
  }

  // ---------------------------------------------------------------------------
  /// Parse and store information about next energy term
  void ParseEnergyTerm(istream &in, Token &token, int nimages, int npsets)
  {
    double weight = Number(in, token);
    if (token == MUL) {
      token = NextToken(in);
    }
    if (token != NAME) {
      cout << endl;
      cerr << "Expected similarity or constraint name after optional weight value" << endl;
      exit(1);
    }

    string             name = token._Name;
    TransformationInfo target_transformation;
    TransformationInfo source_transformation;
    Array<int>         target_index;
    Array<int>         source_index;
    bool               source_is_image = false;

    // Determine possible type of energy term
    //
    // Note: Name may be valid for more than just one type of energy term.
    //       Therefore, do not use "else if" below. Instead, decide which
    //       type of term it is later on when examining the term arguments.
    enum              { SIM,   PDM,   PCM,   CM   , MSDE  };
    bool term_is_[] = { false, false, false, false, false };

    SimilarityMeasure       similarity = _Filter->_SimilarityMeasure;
    PointSetDistanceMeasure pdm        = _Filter->_PointSetDistanceMeasure;
    ConstraintMeasure       constraint;

    // InternalForceTerm enum only available when MIRTK_HAVE_Deformable
    // Hence, use EnergyMeasure enumeration of Common module instead
    EnergyMeasure measure;
    FromString(name.c_str(), measure);
    const string lname = ToLower(name);
    term_is_[PCM ] = (IFT_Begin < measure && measure < IFT_End);
    term_is_[SIM ] = ((lname == "sim") || FromString(name.c_str(), similarity));
    term_is_[PDM ] = ((lname == "pdm") || FromString(name.c_str(), pdm));
    term_is_[CM  ] = (                    FromString(name.c_str(), constraint));
    term_is_[MSDE] = (measure == EM_MeanSquaredDisplacementError);
    if (!term_is_[SIM] && !term_is_[PDM] && !term_is_[PCM] && !term_is_[CM] && !term_is_[MSDE]) {
      cout << endl;
      cerr << "Unknown energy term: " << name << endl;
      exit(1);
    }

    bool default_measure = (term_is_[SIM] && lname == "sim") ||
                           (term_is_[PDM] && lname == "pdm");

    // Custom name
    token = NextToken(in);
    if (token == LB) {
      token = NextToken(in, false);
      name.clear();
      bool var_name = false;
      while (token == NAME  || token == NUMBER || token == SPACE || token == COLON ||
             token == PLUS  || token == MINUS  || token == MUL   || token == DIV   ||
             token == CARET || token == CIRC   || token == LCB   || token == RCB) {
        if      (token == SPACE)  name += ' ';
        else if (token == COLON)  name += ':';
        else if (token == PLUS)   name += '+';
        else if (token == MINUS)  name += '-';
        else if (token == MUL)    name += '*';
        else if (token == DIV)    name += '/';
        else if (token == CARET)  name += '^';
        else if (token == CIRC)   name += 'o';
        else if (token == NUMBER) name += ToString(token._Number);
        else if (token == LCB) {
          if (var_name) {
            cerr << "Expected closing } in custom energy term name" << name << endl;
            exit(1);
          }
          var_name = true;
          name += '{';
        } else if (token == RCB) {
          if (!var_name) {
            cerr << "Found closing } without opening { in custom energy term name: " << name << endl;
            exit(1);
          }
          var_name = false;
          name += '}';
        } else if (token == NAME) {
          name += token._Name;
        } else {
          cerr << "Internal error: Unhandled token in energy term name" << endl;
          exit(1);
        }
        token = NextToken(in, false);
      }
      if (var_name) {
        cerr << "Missing closing } in custom energy term name: " << name << endl;
        exit(1);
      }
      if (token != RB) {
        cout << endl;
        cerr << "Expected closing ] after custom energy term name: " << name << endl;
        exit(1);
      }
      token = NextToken(in);
    }

    // Arguments
    if (token == LP) {
      token = NextToken(in);
      if (token == NAME) {
        // ---------------------------------------------------------------------
        // left composition with transformation as needed for point sets
        if (token._Name == "T") {
          token = NextToken(in);
          if (token == RP) {
            if (!term_is_[CM] && !term_is_[MSDE]) {
              cout << endl;
              cerr << "Transformation identifier cannot be only argument of energy term \"" << name << "\"," << endl;
              cerr << "which is not a (known) transformation constraint (i.e., regularization term)." << endl;
              exit(1);
            }
          } else {
            if (token == CARET) {
              token = NextToken(in);
              if (token != NUMBER && token != MINUS && token != PLUS) {
                cout << endl;
                cerr << "Expected exponent after ^ of transformation composition in argument list of energy term: " << name << endl;
                exit(1);
              }
              target_transformation = Number(in, token);
            } else {
              target_transformation = 1.0;
            }
            if (token != CIRC) {
              cout << endl;
              cerr << "Expected function composition sign (o) after transformation identifier in argument list of energy term: " << name << endl;
              exit(1);
            }
            token = NextToken(in);
            term_is_[CM  ] = false;
            term_is_[MSDE] = false;
          }
        }
        // ---------------------------------------------------------------------
        if ((term_is_[CM] || term_is_[MSDE]) && token == RP) {
          // Term must be a transformation constraint
          term_is_[SIM] = false;
          term_is_[PDM] = false;
          term_is_[PCM] = false;
        // ---------------------------------------------------------------------
        } else if (!(target_index = InputIndex(in, token, name, IMAGE, nimages)).empty()) {
          // Term must be an image similarity
          if (!term_is_[SIM]) {
            cout << endl;
            cerr << "Image identifier cannot be argument of energy term \"" << name << "\"," << endl;
            cerr << "which is not a (known) pairwise image similarity term." << endl;
            exit(1);
          }
          term_is_[PDM]  = false;
          term_is_[CM]   = false;
          term_is_[MSDE] = false;
          // Composition with transformation from left not valid
          if (target_transformation) {
            cout << endl;
            cerr << "Invalid function composition (o) of first image in argument list of fiducial registration error: " << name << endl;
            cerr << "Note that the (source) image is transformed using the notation I2 o T, not T o I2!" << endl;
            exit(1);
          }
          // Target image argument
          if (token == CIRC) {
            token = NextToken(in);
            if (token != NAME || token._Name != "T") {
              cout << endl;
              cerr << "Expected transformation after function composition sign (o) in argument list of similarity term: " << name << endl;
              exit(1);
            }
            token = NextToken(in);
            if (token == CARET) {
              token = NextToken(in);
              if (token != NUMBER && token != MINUS && token != PLUS) {
                cout << endl;
                cerr << "Expected exponent after ^ of transformation composition in argument list of similarity term: " << name << endl;
                exit(1);
              }
              target_transformation = Number(in, token);
            } else {
              target_transformation = 1.0;
            }
          }
          // Source image argument
          if (token != COMMA) {
            cout << endl;
            cerr << "Expected separating , after first image identifier in argument list of similarity term: " << name << endl;
            exit(1);
          }
          token = NextToken(in);
          if (token == NAME) source_index = InputIndex(in, token, name, IMAGE, nimages);
          if (source_index.empty()) {
            cout << endl;
            cerr << "Expected second image identifier after , of argument list of similarity term: " << name << endl;
            exit(1);
          }
          if (token == CIRC) {
            token = NextToken(in);
            if (token != NAME || token._Name != "T") {
              cout << endl;
              cerr << "Expected transformation after function composition sign (o) in argument list of similarity term: " << name << endl;
              exit(1);
            }
            token = NextToken(in);
            if (token == CARET) {
              token = NextToken(in);
              if (token != NUMBER && token != MINUS && token != PLUS) {
                cout << endl;
                cerr << "Expected exponent after ^ of transformation composition in argument list of similarity term: " << name << endl;
                exit(1);
              }
              source_transformation = Number(in, token);
            } else {
              source_transformation = 1.0;
            }
          }
          // By default, transform second image if no explicit transformation was specified
          if (target_transformation.IsIdentity() && source_transformation.IsIdentity()) {
            source_transformation = 1.0;
          }
          source_is_image = true;
        // ---------------------------------------------------------------------
        } else if (!(target_index = InputIndex(in, token, name, POLYDATA, npsets)).empty()) {
          // Term must be a point set distance or constraint
          if (!term_is_[PDM] && !term_is_[PCM]) {
            cout << endl;
            cerr << "Point set identifier cannot be argument of energy term \"" << name << "\"," << endl;
            cerr << "which is neither a (known) pairwise point set distance or internal force term." << endl;
            exit(1);
          }
          // Term cannot be image similarity or transformation constraint
          term_is_[SIM]  = false;
          term_is_[CM]   = false;
          term_is_[MSDE] = false;
          // Target data set argument
          if (token == CIRC) {
            cout << endl;
            cerr << "Invalid function composition (o) of point set in argument list of energy term: " << name << endl;
            cerr << "Note that data sets are transformed using the notation T o S1, not S1 o T!" << endl;
            exit(1);
          }
          if (token == RP) {
            // Term must be a point set constraint (i.e., internal force)
            if (!term_is_[PCM]) {
              cout << endl;
              cerr << "Point set identifier cannot be argument of energy term \"" << name << "\"," << endl;
              cerr << "which is not a (known) internal point set force term." << endl;
              exit(1);
            }
            term_is_[PDM] = false;
            // By default, transform data set if no explicit transformation was specified
            if (target_transformation.IsIdentity()) {
              target_transformation = 1.0;
            }
          } else {
            // Source data set argument
            if (token != COMMA) {
              cout << endl;
              cerr << "Expected separating , after first point set identifier in argument list of point set energy term: " << name << endl;
              exit(1);
            }
            token = NextToken(in);
            if (token != NAME) {
              cout << endl;
              cerr << "Expected transformation or second input identifier after , of argument list of point set energy term: " << name << endl;
              exit(1);
            }
            source_index = InputIndex(in, token, name, POLYDATA, npsets);
            if (source_index.empty()) {
              if (term_is_[PCM]) {
                source_index = InputIndex(in, token, name, IMAGE, nimages);
                if (source_index.empty()) {
                  cout << endl;
                  if (token == NAME && token._Name == "T") {
                    cerr << "Expected second input identifier after , of argument list of point set constraint term: " << name << endl;
                    cerr << "Note that second data set is only used to define reference time frame to which the first" << endl;
                    cerr << "point set is transformed to and therefore no transformation of the second data set is allowed." << endl;
                  } else {
                    cerr << "Expected second input identifier after , of argument list of point set constraint term: " << name << endl;
                  }
                  exit(1);
                }
                source_is_image = true;
              } else {
                if (token != NAME || token._Name != "T") {
                  cout << endl;
                  cerr << "Expected (optionally transformed) second input identifier after , of argument list of point set energy term: " << name << endl;
                  exit(1);
                }
                token = NextToken(in);
                if (token == CARET) {
                  token = NextToken(in);
                  if (token != NUMBER && token != MINUS && token != PLUS) {
                    cout << endl;
                    cerr << "Expected exponent after ^ of transformation composition in argument list of point set energy term: " << name << endl;
                    exit(1);
                  }
                  source_transformation = Number(in, token);
                } else {
                  source_transformation = 1.0;
                }
                if (token != CIRC) {
                  cout << endl;
                  cerr << "Expected function composition sign (o) after transformation identifier in argument list of point set energy term: " << name << endl;
                  exit(1);
                }
                token = NextToken(in);
                if (token == NAME) source_index = InputIndex(in, token, name, POLYDATA, npsets);
                if (source_index.empty()) {
                  cout << endl;
                  cerr << "Expected point set identifier after function composition sign (o) in argument list of point set distance: " << name << endl;
                  exit(1);
                }
              }
            }
            if (token == CIRC) {
              cout << endl;
              cerr << "Invalid function composition (o) of second data set in argument list of point set energy term: " << name << endl;
              if (term_is_[PDM]) {
                cerr << "Note that (source) point sets are transformed using the notation T^-1 o S2, not S2 o T^-1!" << endl;
              }
              exit(1);
            }
            // By default, transform first data set if no explicit transformation was specified
            if (target_transformation.IsIdentity() && source_transformation.IsIdentity()) {
              target_transformation = 1.0;
            }
          }
        // ---------------------------------------------------------------------
        } else {
          cout << endl;
          cerr << "Expected ";
          if (term_is_[SIM]) cerr << "image";
          if (term_is_[PDM] || term_is_[PCM]) {
            if (term_is_[SIM]) cerr << " or ";
            cerr << "point set";
          }
          if (term_is_[CM] || term_is_[MSDE]) {
            if (term_is_[SIM] || term_is_[PDM] || term_is_[PCM]) cerr << " or ";
            cerr << "transformation";
          }
          cerr << " identifier as argument of energy term: " << name << endl;
          exit(1);
        }
      }
      if (token != RP) {
        cout << endl;
        cerr << "Expected closing ) after argument list of energy term: " << name << endl;
        exit(1);
      }
      token = NextToken(in);
    }
    if (target_index.empty() && source_index.empty()) {
      if (!term_is_[CM] && !term_is_[MSDE]) {
        cout << endl;
        cerr << "Expected ";
        if (term_is_[SIM]) cerr << "image";
        if (term_is_[PDM] || term_is_[PCM]) {
          if (term_is_[SIM]) cerr << " or ";
          cerr << "point set";
        }
        cerr << " identifiers as argument of energy term: " << name << endl;
        exit(1);
      }
      term_is_[SIM] = false;
      term_is_[PDM] = false;
      term_is_[PCM] = false;
    }

    // Add energy term structure with parsed information
    if (term_is_[SIM]) {

      if (term_is_[PDM] || term_is_[PCM] || term_is_[CM] || term_is_[MSDE]) {
        cout << endl;
        cerr << "Ambiguous energy term: " << name << endl;
        exit(1);
      }
      if (target_index.size() == 0 || source_index.size() == 0) {
        cout << endl;
        cerr << "Missing image identifier in argument list of similarity term: " << name << endl;
        exit(1);
      }
      if (!source_is_image) {
        cout << endl;
        cerr << "Only images allowed as inputs of the image dissimilarity term: " << name << endl;
        exit(1);
      }

      if (target_index.size() > 1) {
        if (source_index.size() != 1 && source_index.size() != target_index.size()) {
          cout << endl;
          cerr << "Mismatch of image index ranges in argument list of similarity term: " << name << endl;
          exit(1);
        }
        weight /= target_index.size();
        for (size_t t = 0; t < target_index.size(); ++t) {
          size_t s = (source_index.size() == 1 ? 0 : t);
          ImageSimilarityInfo info;
          info._DefaultSign          = default_measure;
          info._Weight               = weight;
          info._Measure              = similarity;
          info._TargetIndex          = target_index[t];
          info._TargetTransformation = target_transformation;
          info._SourceIndex          = source_index[s];
          info._SourceTransformation = source_transformation;
          info._Name                 = TermName(name, info, static_cast<int>(t));
          _Filter->_ImageSimilarityInfo.push_back(info);
        }
      } else {
        if (target_index.size() != 1) {
          cout << endl;
          cerr << "Mismatch of image index ranges in argument list of similarity term: " << name << endl;
          exit(1);
        }
        weight /= source_index.size();
        for (size_t s = 0; s < source_index.size(); ++s) {
          ImageSimilarityInfo info;
          info._DefaultSign          = default_measure;
          info._Weight               = weight;
          info._Measure              = similarity;
          info._TargetIndex          = target_index[0];
          info._TargetTransformation = target_transformation;
          info._SourceIndex          = source_index[s];
          info._SourceTransformation = source_transformation;
          info._Name                 = TermName(name, info, static_cast<int>(s));
          _Filter->_ImageSimilarityInfo.push_back(info);
        }
      }

    }

    if (term_is_[PDM]) {

      if (term_is_[SIM] || term_is_[PCM] || term_is_[CM] || term_is_[MSDE]) {
        cout << endl;
        cerr << "Ambiguous energy term: " << name << endl;
        exit(1);
      }
      if (target_index.size() == 0 || source_index.size() == 0) {
        cout << endl;
        cerr << "Missing point set identifier in argument list of point set distance term: " << name << endl;
        exit(1);
      }
      if (source_is_image) {
        cout << endl;
        cerr << "Only point sets allowed as inputs of the point set distance term: " << name << endl;
        exit(1);
      }

      if (target_index.size() > 1) {
        if (source_index.size() != 1 && source_index.size() != target_index.size()) {
          cout << endl;
          cerr << "Mismatch of point set index ranges in argument list of point set distance term: " << name << endl;
          exit(1);
        }
        weight /= target_index.size();
        for (size_t t = 0; t < target_index.size(); ++t) {
          size_t s = (source_index.size() == 1 ? 0 : t);
          PointSetDistanceInfo info;
          info._DefaultSign          = default_measure;
          info._Weight               = weight;
          info._Measure              = pdm;
          info._TargetIndex          = target_index[t];
          info._TargetTransformation = target_transformation;
          info._SourceIndex          = source_index[s];
          info._SourceTransformation = source_transformation;
          info._Name                 = TermName(name, info, static_cast<int>(t));
          _Filter->_PointSetDistanceInfo.push_back(info);
        }
      } else {
        if (target_index.size() != 1) {
          cout << endl;
          cerr << "Mismatch of point set index ranges in argument list of point set distance term: " << name << endl;
          exit(1);
        }
        weight /= source_index.size();
        for (size_t s = 0; s < source_index.size(); ++s) {
          PointSetDistanceInfo info;
          info._DefaultSign          = default_measure;
          info._Weight               = weight;
          info._Measure              = pdm;
          info._TargetIndex          = target_index[0];
          info._TargetTransformation = target_transformation;
          info._SourceIndex          = source_index[s];
          info._SourceTransformation = source_transformation;
          info._Name                 = TermName(name, info, static_cast<int>(s));
          _Filter->_PointSetDistanceInfo.push_back(info);
        }
      }

    }

    if (term_is_[PCM]) {

      if (term_is_[SIM] || term_is_[PDM] || term_is_[CM] || term_is_[MSDE]) {
        cout << endl;
        cerr << "Ambiguous energy term: " << name << endl;
        exit(1);
      }
      if (target_index.size() == 0) {
        cout << endl;
        cerr << "Missing point set identifier in argument list of internal point set force term: " << name << endl;
        exit(1);
      }

      if (source_index.size() > 0) {
        if (target_index.size() > 1) {
          if (source_index.size() != 1 && source_index.size() != target_index.size()) {
            cout << endl;
            cerr << "Mismatch of input index ranges in argument list of point set constraint term: " << name << endl;
            exit(1);
          }
          weight /= target_index.size();
          for (size_t t = 0; t < target_index.size(); ++t) {
            size_t s = (source_index.size() == 1 ? 0 : t);
            PointSetConstraintInfo info;
            info._Weight               = weight;
            info._Measure              = measure;
            info._PointSetIndex        = target_index[t];
            info._Transformation       = target_transformation;
            if (source_is_image) {
              info._RefPointSetIndex   = -1;
              info._RefImageIndex      = source_index[s];
            } else {
              info._RefPointSetIndex   = source_index[s];
              info._RefImageIndex      = -1;
            }
            info._Name = TermName(name, info, static_cast<int>(t));
            _Filter->_PointSetConstraintInfo.push_back(info);
          }
        } else {
          if (target_index.size() != 1) {
            cout << endl;
            cerr << "Mismatch of input index ranges in argument list of point set constraint term: " << name << endl;
            exit(1);
          }
          weight /= source_index.size();
          for (size_t s = 0; s < source_index.size(); ++s) {
            PointSetConstraintInfo info;
            info._Weight               = weight;
            info._Measure              = measure;
            info._PointSetIndex        = target_index[0];
            info._Transformation       = target_transformation;
            if (source_is_image) {
              info._RefPointSetIndex   = -1;
              info._RefImageIndex      = source_index[s];
            } else {
              info._RefPointSetIndex   = source_index[s];
              info._RefImageIndex      = -1;
            }
            info._Name = TermName(name, info, static_cast<int>(s));
            _Filter->_PointSetConstraintInfo.push_back(info);
          }
        }
      } else {
        weight /= target_index.size();
        for (size_t t = 0; t < target_index.size(); ++t) {
          PointSetConstraintInfo info;
          info._Weight           = weight;
          info._Measure          = measure;
          info._PointSetIndex    = target_index[t];
          info._Transformation   = target_transformation;
          info._RefPointSetIndex = -1;
          info._RefImageIndex    = -1;
          // {i} substituted by GenericRegistrationFilter::AddPointSetConstraint
          info._Name = TermName(name, info, -1);
          _Filter->_PointSetConstraintInfo.push_back(info);
        }
      }

    }

    if (term_is_[CM]) {

      if (term_is_[SIM] || term_is_[PDM] || term_is_[PCM] || term_is_[MSDE]) {
        cout << endl;
        cerr << "Ambiguous energy term: " << name << endl;
        exit(1);
      }

      ConstraintInfo info;
      info._Weight  = weight;
      info._Measure = constraint;
      info._Name    = name;
      _Filter->_ConstraintInfo.push_back(info);

    }

    if (term_is_[MSDE]) {

      if (term_is_[SIM] || term_is_[PDM] || term_is_[PCM] || term_is_[CM]) {
        cout << endl;
        cerr << "Ambiguous energy term: " << name << endl;
        exit(1);
      }

      _Filter->TargetTransformationErrorName(name);
      _Filter->TargetTransformationErrorWeight(weight);

    }

    if (!term_is_[SIM] && !term_is_[PDM] && !term_is_[PCM] && !term_is_[CM] && !term_is_[MSDE]) {
      cout << endl;
      cerr << "Unknown enery term: " << name << endl;
      exit(1);
    }
  }

public:

  // ---------------------------------------------------------------------------
  /// Constructor
  RegistrationEnergyParser(GenericRegistrationFilter *filter)
  :
    _Filter(filter)
  {}

  // ---------------------------------------------------------------------------
  /// Parse energy function
  void ParseEnergyFormula(const string &energy_formula, int nimages = -1, int npsets = -1)
  {
    istringstream formula(energy_formula);
    Token              token  (END);

    _Filter->_ImageSimilarityInfo   .clear();
    _Filter->_PointSetDistanceInfo  .clear();
    _Filter->_PointSetConstraintInfo.clear();
    _Filter->_ConstraintInfo        .clear();

    if (nimages < 0) nimages = _Filter->NumberOfImages();
    if (npsets  < 0) npsets  = _Filter->NumberOfPointSets();

    token = NextToken(formula);
    while (token != END) {
      ParseEnergyTerm(formula, token, nimages, npsets);
    }
  }

}; // RegistrationEnergyParser


} // namespace mirtk
