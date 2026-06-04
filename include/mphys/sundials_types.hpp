#pragma once

#include <utility>

#include <ida/ida.h>
#include <kinsol/kinsol.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_band.h>

namespace mphys {

struct SunContext {
  SUNContext ctx = nullptr;
  SunContext() { SUNContext_Create(SUN_COMM_NULL, &ctx); }
  ~SunContext() {
    if (ctx) SUNContext_Free(&ctx);
  }
  SunContext(const SunContext&) = delete;
  SunContext& operator=(const SunContext&) = delete;
  SunContext(SunContext&& other) noexcept : ctx(std::exchange(other.ctx, nullptr)) {}
  SunContext& operator=(SunContext&& other) noexcept {
    if (this != &other) {
      if (ctx) SUNContext_Free(&ctx);
      ctx = std::exchange(other.ctx, nullptr);
    }
    return *this;
  }
  operator SUNContext() const { return ctx; }
};

struct SunVector {
  N_Vector v = nullptr;
  SunVector() = default;
  SunVector(sunindextype n, SUNContext ctx) : v(N_VNew_Serial(n, ctx)) {}
  ~SunVector() {
    if (v) N_VDestroy(v);
  }
  SunVector(const SunVector&) = delete;
  SunVector& operator=(const SunVector&) = delete;
  SunVector(SunVector&& other) noexcept : v(std::exchange(other.v, nullptr)) {}
  SunVector& operator=(SunVector&& other) noexcept {
    if (this != &other) {
      if (v) N_VDestroy(v);
      v = std::exchange(other.v, nullptr);
    }
    return *this;
  }
  operator N_Vector() const { return v; }
  double* Data() { return N_VGetArrayPointer(v); }
  const double* Data() const { return N_VGetArrayPointer(v); }
};

struct SunMatrix {
  SUNMatrix m = nullptr;
  SunMatrix() = default;
  SunMatrix(sunindextype n, sunindextype mu, sunindextype ml, SUNContext ctx)
      : m(SUNBandMatrix(n, mu, ml, ctx)) {}
  ~SunMatrix() {
    if (m) SUNMatDestroy(m);
  }
  SunMatrix(const SunMatrix&) = delete;
  SunMatrix& operator=(const SunMatrix&) = delete;
  SunMatrix(SunMatrix&& other) noexcept : m(std::exchange(other.m, nullptr)) {}
  SunMatrix& operator=(SunMatrix&& other) noexcept {
    if (this != &other) {
      if (m) SUNMatDestroy(m);
      m = std::exchange(other.m, nullptr);
    }
    return *this;
  }
  operator SUNMatrix() const { return m; }
};

struct SunLinearSolver {
  SUNLinearSolver ls = nullptr;
  SunLinearSolver() = default;
  SunLinearSolver(N_Vector y, SUNMatrix a, SUNContext ctx)
      : ls(SUNLinSol_Band(y, a, ctx)) {}
  // Matrix-free Krylov solver (GMRES). pretype is a SUN_PREC_* constant; maxl
  // is the maximum Krylov subspace dimension.
  SunLinearSolver(N_Vector y, int pretype, int maxl, SUNContext ctx)
      : ls(SUNLinSol_SPGMR(y, pretype, maxl, ctx)) {}
  ~SunLinearSolver() {
    if (ls) SUNLinSolFree(ls);
  }
  SunLinearSolver(const SunLinearSolver&) = delete;
  SunLinearSolver& operator=(const SunLinearSolver&) = delete;
  SunLinearSolver(SunLinearSolver&& other) noexcept
      : ls(std::exchange(other.ls, nullptr)) {}
  SunLinearSolver& operator=(SunLinearSolver&& other) noexcept {
    if (this != &other) {
      if (ls) SUNLinSolFree(ls);
      ls = std::exchange(other.ls, nullptr);
    }
    return *this;
  }
  operator SUNLinearSolver() const { return ls; }
};

struct IdaMem {
  void* mem = nullptr;
  IdaMem() = default;
  explicit IdaMem(SUNContext ctx) : mem(IDACreate(ctx)) {}
  ~IdaMem() {
    if (mem) IDAFree(&mem);
  }
  IdaMem(const IdaMem&) = delete;
  IdaMem& operator=(const IdaMem&) = delete;
  IdaMem(IdaMem&& other) noexcept : mem(std::exchange(other.mem, nullptr)) {}
  IdaMem& operator=(IdaMem&& other) noexcept {
    if (this != &other) {
      if (mem) IDAFree(&mem);
      mem = std::exchange(other.mem, nullptr);
    }
    return *this;
  }
  operator void*() const { return mem; }
};

struct KinMem {
  void* mem = nullptr;
  KinMem() = default;
  explicit KinMem(SUNContext ctx) : mem(KINCreate(ctx)) {}
  ~KinMem() {
    if (mem) KINFree(&mem);
  }
  KinMem(const KinMem&) = delete;
  KinMem& operator=(const KinMem&) = delete;
  KinMem(KinMem&& other) noexcept : mem(std::exchange(other.mem, nullptr)) {}
  KinMem& operator=(KinMem&& other) noexcept {
    if (this != &other) {
      if (mem) KINFree(&mem);
      mem = std::exchange(other.mem, nullptr);
    }
    return *this;
  }
  operator void*() const { return mem; }
};

}  // namespace mphys
