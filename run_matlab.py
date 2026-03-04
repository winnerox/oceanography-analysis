import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import subprocess
import threading
import os
import sys
from pathlib import Path
import queue
import time
import re

class ScriptRunnerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("MATLAB脚本批量运行工具")
        self.root.geometry("1000x750")
        
        # 设置样式
        self.root.configure(bg='#f0f0f0')
        
        # 存储脚本路径的列表
        self.script_paths = []
        
        # 存储运行中的进程
        self.processes = []
        
        # 输出队列
        self.output_queue = queue.Queue()
        
        # 创建GUI组件
        self.create_widgets()
        
        # 开始处理输出队列
        self.process_output_queue()
        
        # 设置拖拽功能
        self.setup_drag_drop()
        
    def setup_drag_drop(self):
        """设置拖拽功能"""
        # 这需要使用tkinterdnd2库，但为了简化，我们用按钮模拟
        # 实际拖拽功能需要安装：pip install tkinterdnd2
        pass
        
    def create_widgets(self):
        # 标题
        title_label = tk.Label(self.root, text="MATLAB脚本批量运行工具", 
                               font=("Arial", 16, "bold"), bg='#f0f0f0', fg='#333')
        title_label.pack(pady=10)
        
        # 拖拽上传区域
        self.create_drop_area()
        
        # 主内容区域（左右分栏）
        main_frame = tk.Frame(self.root, bg='#f0f0f0')
        main_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=10)
        
        # 左侧：脚本列表
        left_frame = tk.Frame(main_frame, bg='#f0f0f0', width=400)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 10))
        left_frame.pack_propagate(False)
        
        self.create_script_list(left_frame)
        
        # 右侧：运行设置
        right_frame = tk.Frame(main_frame, bg='#f0f0f0', width=300)
        right_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=(10, 0))
        right_frame.pack_propagate(False)
        
        self.create_settings_panel(right_frame)
        
        # 控制面板
        self.create_control_panel()
        
        # 输出区域
        self.create_output_area()
        
        # 状态栏
        self.status_bar = tk.Label(self.root, text="就绪", bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
    def create_drop_area(self):
        drop_frame = tk.Frame(self.root, bg='#e0e0e0', height=80, relief=tk.GROOVE, bd=2)
        drop_frame.pack(fill=tk.X, padx=20, pady=10)
        drop_frame.pack_propagate(False)
        
        # 拖拽提示
        self.drop_label = tk.Label(drop_frame, 
                                   text="📁 点击此处选择MATLAB脚本文件\n或直接将文件拖拽到此窗口",
                                   font=("Arial", 11),
                                   bg='#e0e0e0', fg='#666',
                                   cursor="hand2")
        self.drop_label.pack(expand=True, fill=tk.BOTH)
        
        # 绑定点击事件
        self.drop_label.bind('<Button-1>', lambda e: self.select_files())
        
        # 说明文本
        info_text = tk.Label(self.root, text="支持 .m 文件，可多选", 
                            font=("Arial", 9), bg='#f0f0f0', fg='#999')
        info_text.pack()
        
    def create_script_list(self, parent):
        tk.Label(parent, text="📋 已选择脚本列表", font=("Arial", 11, "bold"), 
                bg='#f0f0f0').pack(anchor=tk.W, pady=(0, 5))
        
        # 创建列表和滚动条
        listbox_frame = tk.Frame(parent, bg='#f0f0f0')
        listbox_frame.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = tk.Scrollbar(listbox_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.script_listbox = tk.Listbox(listbox_frame, 
                                         selectmode=tk.EXTENDED,
                                         yscrollcommand=scrollbar.set,
                                         height=12,
                                         font=("Consolas", 10),
                                         bg='white',
                                         selectbackground='#4CAF50')
        self.script_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        scrollbar.config(command=self.script_listbox.yview)
        
        # 列表下方的按钮
        btn_frame = tk.Frame(parent, bg='#f0f0f0')
        btn_frame.pack(fill=tk.X, pady=10)
        
        tk.Button(btn_frame, text="➕ 添加脚本", command=self.select_files,
                 bg='#4CAF50', fg='white', font=("Arial", 9),
                 width=10).pack(side=tk.LEFT, padx=2)
        tk.Button(btn_frame, text="➖ 移除选中", command=self.remove_selected,
                 bg='#f44336', fg='white', font=("Arial", 9),
                 width=10).pack(side=tk.LEFT, padx=2)
        tk.Button(btn_frame, text="🗑️ 清空列表", command=self.clear_list,
                 bg='#ff9800', fg='white', font=("Arial", 9),
                 width=10).pack(side=tk.LEFT, padx=2)
        
    def create_settings_panel(self, parent):
        # 运行模式设置
        tk.Label(parent, text="⚙️ 运行设置", font=("Arial", 11, "bold"),
                bg='#f0f0f0').pack(anchor=tk.W, pady=(0, 10))
        
        # 模式选择框架
        mode_frame = tk.Frame(parent, bg='#ffffff', relief=tk.GROOVE, bd=1)
        mode_frame.pack(fill=tk.X, pady=5)
        
        self.run_mode = tk.StringVar(value="sequential")
        
        modes = [
            ("🔂 按顺序运行所有", "sequential"),
            ("🔄 两两并行运行", "parallel_two"),
            ("⚡ 全部并行运行", "parallel_all")
        ]
        
        for text, mode in modes:
            rb = tk.Radiobutton(mode_frame, text=text, variable=self.run_mode, 
                               value=mode, bg='white', font=("Arial", 9),
                               anchor=tk.W, padx=10)
            rb.pack(fill=tk.X, pady=5)
        
        # 分隔线
        ttk.Separator(parent, orient='horizontal').pack(fill=tk.X, pady=15)
        
        # MATLAB路径设置
        tk.Label(parent, text="📌 MATLAB路径:", font=("Arial", 10, "bold"),
                bg='#f0f0f0').pack(anchor=tk.W, pady=(0, 5))
        
        path_frame = tk.Frame(parent, bg='#f0f0f0')
        path_frame.pack(fill=tk.X, pady=2)
        
        self.matlab_path = tk.StringVar(value="matlab")
        path_entry = tk.Entry(path_frame, textvariable=self.matlab_path, 
                             font=("Consolas", 9), bg='white')
        path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))
        
        tk.Button(path_frame, text="浏览", command=self.browse_matlab,
                 bg='#607d8b', fg='white', font=("Arial", 8),
                 width=6).pack(side=tk.RIGHT)
        
        # 帮助信息
        help_text = "💡 提示：如果matlab命令已在系统PATH中，\n   直接输入'matlab'即可"
        tk.Label(parent, text=help_text, font=("Arial", 8),
                bg='#f0f0f0', fg='#666', justify=tk.LEFT).pack(anchor=tk.W, pady=10)
        
        # 额外选项
        tk.Label(parent, text="🔧 高级选项", font=("Arial", 10, "bold"),
                bg='#f0f0f0').pack(anchor=tk.W, pady=(10, 5))
        
        self.show_output_var = tk.BooleanVar(value=True)
        tk.Checkbutton(parent, text="显示实时输出", variable=self.show_output_var,
                      bg='#f0f0f0', font=("Arial", 9)).pack(anchor=tk.W)
        
    def create_control_panel(self):
        control_frame = tk.Frame(self.root, bg='#e0e0e0', height=60)
        control_frame.pack(fill=tk.X, padx=20, pady=10)
        control_frame.pack_propagate(False)
        
        # 运行按钮
        self.run_btn = tk.Button(control_frame, text="▶ 运行选中脚本", 
                                 command=self.run_scripts,
                                 bg='#2196F3', fg='white', font=("Arial", 12, "bold"),
                                 width=15, height=1)
        self.run_btn.pack(side=tk.LEFT, padx=20)
        
        # 停止按钮
        self.stop_btn = tk.Button(control_frame, text="⏹ 停止运行", 
                                  command=self.stop_scripts,
                                  bg='#9e9e9e', fg='white', font=("Arial", 12),
                                  width=15, height=1, state=tk.DISABLED)
        self.stop_btn.pack(side=tk.LEFT, padx=10)
        
        # 清空输出按钮
        tk.Button(control_frame, text="🗑️ 清空输出", 
                 command=self.clear_output,
                 bg='#607d8b', fg='white', font=("Arial", 10),
                 width=10).pack(side=tk.RIGHT, padx=20)
        
    def create_output_area(self):
        output_frame = tk.Frame(self.root)
        output_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=(0, 10))
        
        # 标题和计数器
        header_frame = tk.Frame(output_frame)
        header_frame.pack(fill=tk.X)
        
        tk.Label(header_frame, text="📊 运行输出", font=("Arial", 11, "bold")).pack(side=tk.LEFT)
        
        self.output_counter = tk.Label(header_frame, text="", font=("Arial", 9), fg='#666')
        self.output_counter.pack(side=tk.RIGHT)
        
        # 创建文本框和滚动条
        text_frame = tk.Frame(output_frame)
        text_frame.pack(fill=tk.BOTH, expand=True)
        
        scrollbar = tk.Scrollbar(text_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.output_text = tk.Text(text_frame, 
                                    wrap=tk.WORD,
                                    yscrollcommand=scrollbar.set,
                                    height=15,
                                    font=("Consolas", 10),
                                    bg='#1e1e1e', fg='#d4d4d4')
        self.output_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        scrollbar.config(command=self.output_text.yview)
        
        # 配置文本标签颜色
        self.output_text.tag_config("info", foreground="#4ec9b0")
        self.output_text.tag_config("error", foreground="#f48771")
        self.output_text.tag_config("success", foreground="#6a9955")
        self.output_text.tag_config("header", foreground="#569cd6", font=("Consolas", 11, "bold"))
        self.output_text.tag_config("stdout", foreground="#d4d4d4")
        self.output_text.tag_config("stderr", foreground="#f48771")
        
    def select_files(self):
        files = filedialog.askopenfilenames(
            title="选择MATLAB脚本文件",
            filetypes=[("MATLAB文件", "*.m"), ("所有文件", "*.*")]
        )
        
        added = 0
        for file in files:
            if file not in self.script_paths:
                self.script_paths.append(file)
                self.script_listbox.insert(tk.END, os.path.basename(file))
                added += 1
        
        self.update_status(f"已添加 {added} 个脚本，当前共 {len(self.script_paths)} 个")
        
    def remove_selected(self):
        selected = self.script_listbox.curselection()
        removed = len(selected)
        for i in reversed(selected):
            self.script_listbox.delete(i)
            del self.script_paths[i]
        
        self.update_status(f"已移除 {removed} 个，当前共 {len(self.script_paths)} 个")
        
    def clear_list(self):
        self.script_listbox.delete(0, tk.END)
        self.script_paths.clear()
        self.update_status("已清空脚本列表")
        
    def browse_matlab(self):
        path = filedialog.askopenfilename(
            title="选择MATLAB可执行文件",
            filetypes=[("可执行文件", "matlab.exe"), ("所有文件", "*.*")]
        )
        if path:
            # 转换为Windows路径格式
            path = path.replace('/', '\\')
            self.matlab_path.set(path)
            
    def update_status(self, message):
        self.status_bar.config(text=f"📌 {message}")
        self.root.update()
        
    def log_output(self, message, tag="info"):
        """在输出区域添加日志"""
        self.output_queue.put((message, tag))
        
    def process_output_queue(self):
        """处理输出队列"""
        try:
            while True:
                message, tag = self.output_queue.get_nowait()
                self.output_text.insert(tk.END, message + "\n", tag)
                self.output_text.see(tk.END)
                self.output_text.update_idletasks()
                
                # 更新计数器
                line_count = int(self.output_text.index('end-1c').split('.')[0])
                self.output_counter.config(text=f"共 {line_count} 行")
        except queue.Empty:
            pass
        finally:
            self.root.after(100, self.process_output_queue)
        
    def clear_output(self):
        self.output_text.delete(1.0, tk.END)
        self.output_counter.config(text="")
        
    def run_scripts(self):
        if not self.script_paths:
            messagebox.showwarning("警告", "请先添加要运行的脚本")
            return
        
        # 禁用运行按钮，启用停止按钮
        self.run_btn.config(state=tk.DISABLED)
        self.stop_btn.config(state=tk.NORMAL, bg='#f44336')
        
        # 清空旧输出
        self.clear_output()
        
        # 在新线程中运行脚本
        thread = threading.Thread(target=self._run_scripts_thread)
        thread.daemon = True
        thread.start()
        
    def _run_scripts_thread(self):
        """在后台线程中运行脚本"""
        mode = self.run_mode.get()
        matlab_cmd = self.matlab_path.get()
        
        self.log_output("="*70, "header")
        self.log_output(f"🚀 开始运行脚本 (模式: {mode})", "header")
        self.log_output(f"📁 脚本数量: {len(self.script_paths)}", "info")
        self.log_output("="*70, "header")
        
        # 验证MATLAB命令
        if not self._verify_matlab(matlab_cmd):
            self.root.after(0, self._run_completed)
            return
        
        try:
            if mode == "sequential":
                self._run_sequential(matlab_cmd)
            elif mode == "parallel_two":
                self._run_parallel_two(matlab_cmd)
            else:  # parallel_all
                self._run_parallel_all(matlab_cmd)
        except Exception as e:
            self.log_output(f"❌ 运行出错: {str(e)}", "error")
        
        # 运行完成
        self.root.after(0, self._run_completed)
        
    def _verify_matlab(self, matlab_cmd):
        """验证MATLAB命令是否可用"""
        try:
            # 测试MATLAB命令
            cmd = [matlab_cmd, '-batch', "disp('MATLAB is working');"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                self.log_output("✅ MATLAB命令验证成功", "success")
                return True
            else:
                self.log_output(f"❌ MATLAB命令验证失败: {result.stderr}", "error")
                return False
        except Exception as e:
            self.log_output(f"❌ MATLAB命令不可用: {str(e)}", "error")
            return False
        
    def _run_sequential(self, matlab_cmd):
        """按顺序运行"""
        for i, script in enumerate(self.script_paths):
            self.log_output(f"\n{'='*50}", "header")
            self.log_output(f"📌 [{i+1}/{len(self.script_paths)}] 运行: {os.path.basename(script)}", "header")
            self.log_output(f"{'='*50}", "header")
            self._run_single_script(script, matlab_cmd)
            
    def _run_parallel_two(self, matlab_cmd):
        """两两并行运行"""
        import time
        
        for i in range(0, len(self.script_paths), 2):
            batch = self.script_paths[i:i+2]
            processes = []
            
            self.log_output(f"\n{'='*50}", "header")
            self.log_output(f"📌 运行批次 {i//2 + 1}: {[os.path.basename(s) for s in batch]}", "header")
            self.log_output(f"{'='*50}", "header")
            
            # 启动当前批次的所有脚本
            for script in batch:
                self.log_output(f"🚀 启动: {os.path.basename(script)}", "info")
                process = self._run_script_process(script, matlab_cmd)
                if process:
                    processes.append(process)
                    self.processes.append(process)
            
            # 等待当前批次的所有进程完成
            for process in processes:
                process.wait()
                
    def _run_parallel_all(self, matlab_cmd):
        """全部并行运行"""
        processes = []
        
        self.log_output(f"\n{'='*50}", "header")
        self.log_output(f"📌 启动所有脚本 ({len(self.script_paths)}个)", "header")
        self.log_output(f"{'='*50}", "header")
        
        # 启动所有脚本
        for script in self.script_paths:
            self.log_output(f"🚀 启动: {os.path.basename(script)}", "info")
            process = self._run_script_process(script, matlab_cmd)
            if process:
                processes.append(process)
                self.processes.append(process)
        
        # 等待所有进程完成
        for process in processes:
            process.wait()
            
    def _run_script_process(self, script_path, matlab_cmd):
        """创建并返回一个运行脚本的进程"""
        try:
            # 构建命令
            script_path = script_path.replace('/', '\\')
            cmd = [matlab_cmd, '-batch', f"try, run('{script_path}'), catch ME, disp(ME.message), end"]
            
            # 创建进程
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.PIPE,
                text=True,
                bufsize=1,
                universal_newlines=True,
                creationflags=subprocess.CREATE_NO_WINDOW if sys.platform == 'win32' else 0
            )
            
            # 启动线程读取输出
            threading.Thread(target=self._read_output, args=(process, os.path.basename(script_path)), daemon=True).start()
            
            return process
            
        except Exception as e:
            self.log_output(f"❌ 启动失败 {os.path.basename(script_path)}: {str(e)}", "error")
            return None
                
    def _run_single_script(self, script_path, matlab_cmd):
        """运行单个脚本并实时显示输出"""
        process = self._run_script_process(script_path, matlab_cmd)
        if process:
            process.wait()
            
    def _read_output(self, process, script_name):
        """读取进程输出"""
        try:
            # 读取stdout
            for line in iter(process.stdout.readline, ''):
                if line:
                    line = line.strip()
                    if line:
                        self.log_output(f"[{script_name}] {line}", "stdout")
            
            # 读取stderr
            for line in iter(process.stderr.readline, ''):
                if line:
                    line = line.strip()
                    if line:
                        self.log_output(f"[{script_name}] ⚠️ {line}", "stderr")
            
            # 检查返回码
            returncode = process.wait()
            if returncode == 0:
                self.log_output(f"✅ [{script_name}] 运行完成", "success")
            else:
                self.log_output(f"❌ [{script_name}] 运行失败 (返回码: {returncode})", "error")
                
        except Exception as e:
            self.log_output(f"❌ [{script_name}] 读取输出时出错: {str(e)}", "error")
            
    def _run_completed(self):
        """运行完成后的处理"""
        self.processes.clear()
        self.run_btn.config(state=tk.NORMAL)
        self.stop_btn.config(state=tk.DISABLED, bg='#9e9e9e')
        self.log_output("\n" + "="*70, "header")
        self.log_output("🏁 所有脚本运行完成", "header")
        self.log_output("="*70, "header")
        self.update_status("运行完成")
        
    def stop_scripts(self):
        """停止所有运行中的脚本"""
        for process in self.processes:
            try:
                process.terminate()
            except:
                pass
        self.processes.clear()
        self.log_output("\n⚠️ 用户终止运行", "error")
        self._run_completed()

def main():
    root = tk.Tk()
    app = ScriptRunnerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()