#ifndef MODULE_REGISTRATION_H
#define MODULE_REGISTRATION_H

class Calculator;

/**
 * @brief 注册标准数学和系统模块
 * @param calculator 计算器实例指针
 */
void register_standard_modules(Calculator* calculator);

#endif
