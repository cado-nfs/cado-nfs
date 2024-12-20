#ifndef TABLE_F_PY_PHI__F_DEG4_S02_C4V4_H1_PY_S20_F01_H_
#define TABLE_F_PY_PHI__F_DEG4_S02_C4V4_H1_PY_S20_F01_H_

#include "gfpkdlpolyselect.h"

// pragma no prototypes

/* table of polynomials f preselected with good properties. */

fPyphi_t table_f_Py_phi__f_deg4_s02_C4V4_h1_Py_s20_f01[8] = { \
  { /* f= x^4 - x^3 + x^2 - x + 1, Py= x^2 + x - 1, varphi= -y*X + X^2 - X + 1 */
    {/*f     */ 1, -1, 1, -1, 1, }, \
    {/*Py    */ -1, 1, 1, }, \
    {/*varphi*/ { 1, 0, }, { -1, -1, }, { 1, 0, }, }, \
  }, \
  { /* f= x^4 + x^3 + x^2 + x + 1, Py= x^2 - x - 1, varphi= -y*X + X^2 + X + 1 */
    {/*f     */ 1, 1, 1, 1, 1, }, \
    {/*Py    */ -1, -1, 1, }, \
    {/*varphi*/ { 1, 0, }, { 1, -1, }, { 1, 0, }, }, \
  }, \
  { /* f= x^4 - x^2 + 1, Py= x^2 - 3, varphi= -y*X + X^2 + 1 */
    {/*f     */ 1, 0, -1, 0, 1, }, \
    {/*Py    */ -3, 0, 1, }, \
    {/*varphi*/ { 1, 0, }, { 0, -1, }, { 1, 0, }, }, \
  }, \
  { /* f= x^4 - x^3 + 2*x^2 + x + 1, Py= x^2 - 3*x + 1, varphi= -y*X + y + X^2 + X */
    {/*f     */ 1, 1, 2, -1, 1, }, \
    {/*Py    */ 1, -3, 1, }, \
    {/*varphi*/ { 0, 1, }, { 1, -1, }, { 1, 0, }, }, \
  }, \
  { /* f= x^4 + x^3 + 2*x^2 - x + 1, Py= x^2 - 3*x + 1, varphi= -y*X - y + X^2 + 2*X + 3 */
    {/*f     */ 1, -1, 2, 1, 1, }, \
    {/*Py    */ 1, -3, 1, }, \
    {/*varphi*/ { 3, -1, }, { 2, -1, }, { 1, 0, }, }, \
  }, \
  { /* f= x^4 + 1, Py= x^2 - 2, varphi= -y*X + X^2 + 1 */
    {/*f     */ 1, 0, 0, 0, 1, }, \
    {/*Py    */ -2, 0, 1, }, \
    {/*varphi*/ { 1, 0, }, { 0, -1, }, { 1, 0, }, }, \
  }, \
  { /* f= x^4 + 3*x^2 + 1, Py= x^2 - 3*x + 1, varphi= -y + X^2 + 3 */
    {/*f     */ 1, 0, 3, 0, 1, }, \
    {/*Py    */ 1, -3, 1, }, \
    {/*varphi*/ { 3, -1, }, { 0, 0, }, { 1, 0, }, }, \
  }, \
  { /* f= x^4 + 9*x^2 + 1, Py= x^2 - 9*x + 1, varphi= -y + X^2 + 9 */
    {/*f     */ 1, 0, 9, 0, 1, }, \
    {/*Py    */ 1, -9, 1, }, \
    {/*varphi*/ { 9, -1, }, { 0, 0, }, { 1, 0, }, }, \
  }, \
}; 


const fPyphi_poly_t fPyphi_4 = {
  4, 
  2,
  2,
  8,
  table_f_Py_phi__f_deg4_s02_C4V4_h1_Py_s20_f01,
};

#endif
