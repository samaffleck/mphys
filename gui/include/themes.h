#pragma once

#include "imgui.h"

namespace Themes {

  inline void SetDarkTheme() {
    ImGuiStyle& style = ImGui::GetStyle();
    ImVec4* colors = style.Colors;

    // Primary background
    colors[ImGuiCol_WindowBg] = ImVec4(0.07f, 0.07f, 0.09f, 1.00f);  // #131318
    colors[ImGuiCol_MenuBarBg] = ImVec4(0.10f, 0.10f, 0.12f, 1.00f); // Slightly lighter than window bg for contrast
    colors[ImGuiCol_ChildBg] = ImVec4(0.10f, 0.10f, 0.12f, 1.00f); 
    colors[ImGuiCol_PopupBg] = ImVec4(0.18f, 0.18f, 0.22f, 1.00f);
    
    // Progress bar
    colors[ImGuiCol_PlotHistogram] = ImVec4(0.55f, 0.70f, 0.95f, 1.00f);

    // Headers
    colors[ImGuiCol_Header] = ImVec4(0.18f, 0.18f, 0.22f, 1.00f);
    colors[ImGuiCol_HeaderHovered] = ImVec4(0.30f, 0.30f, 0.40f, 1.00f);
    colors[ImGuiCol_HeaderActive] = ImVec4(0.25f, 0.25f, 0.35f, 1.00f);

    // Buttons
    colors[ImGuiCol_Button] = ImVec4(0.20f, 0.22f, 0.27f, 1.00f);
    colors[ImGuiCol_ButtonHovered] = ImVec4(0.30f, 0.32f, 0.40f, 1.00f);
    colors[ImGuiCol_ButtonActive] = ImVec4(0.35f, 0.38f, 0.50f, 1.00f);

    // Frame BG
    colors[ImGuiCol_FrameBg] = ImVec4(0.15f, 0.15f, 0.18f, 1.00f);
    colors[ImGuiCol_FrameBgHovered] = ImVec4(0.22f, 0.22f, 0.27f, 1.00f);
    colors[ImGuiCol_FrameBgActive] = ImVec4(0.25f, 0.25f, 0.30f, 1.00f);

    // Tabs
    colors[ImGuiCol_Tab] = ImVec4(0.18f, 0.18f, 0.22f, 1.00f);
    colors[ImGuiCol_TabHovered] = ImVec4(0.35f, 0.35f, 0.50f, 1.00f);
    colors[ImGuiCol_TabActive] = ImVec4(0.25f, 0.25f, 0.38f, 1.00f);
    colors[ImGuiCol_TabUnfocused] = ImVec4(0.13f, 0.13f, 0.17f, 1.00f);
    colors[ImGuiCol_TabUnfocusedActive] = ImVec4(0.20f, 0.20f, 0.25f, 1.00f);

    // Title
    colors[ImGuiCol_TitleBg] = ImVec4(0.12f, 0.12f, 0.15f, 1.00f);
    colors[ImGuiCol_TitleBgActive] = ImVec4(0.15f, 0.15f, 0.20f, 1.00f);
    colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.10f, 0.10f, 0.12f, 1.00f);

    // Borders
    colors[ImGuiCol_Border] = ImVec4(0.20f, 0.20f, 0.25f, 0.50f);
    colors[ImGuiCol_BorderShadow] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);

    // Text
    colors[ImGuiCol_Text] = ImVec4(0.90f, 0.90f, 0.95f, 1.00f);
    colors[ImGuiCol_TextDisabled] = ImVec4(0.50f, 0.50f, 0.55f, 1.00f);
    colors[ImGuiCol_InputTextCursor] = ImVec4(0.90f, 0.90f, 0.95f, 1.00f); // Match text color

    // Highlights
    colors[ImGuiCol_CheckMark] = ImVec4(0.50f, 0.70f, 1.00f, 1.00f);
    colors[ImGuiCol_SliderGrab] = ImVec4(0.50f, 0.70f, 1.00f, 1.00f);
    colors[ImGuiCol_SliderGrabActive] = ImVec4(0.60f, 0.80f, 1.00f, 1.00f);
    colors[ImGuiCol_ResizeGrip] = ImVec4(0.50f, 0.70f, 1.00f, 0.50f);
    colors[ImGuiCol_ResizeGripHovered] = ImVec4(0.60f, 0.80f, 1.00f, 0.75f);
    colors[ImGuiCol_ResizeGripActive] = ImVec4(0.70f, 0.90f, 1.00f, 1.00f);

    // Scrollbar
    colors[ImGuiCol_ScrollbarBg] = ImVec4(0.10f, 0.10f, 0.12f, 1.00f);
    colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.30f, 0.30f, 0.35f, 1.00f);
    colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.40f, 0.40f, 0.50f, 1.00f);
    colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.45f, 0.45f, 0.55f, 1.00f);

    // Table colors for dark theme
    colors[ImGuiCol_TableHeaderBg]     = ImVec4(0.12f, 0.12f, 0.15f, 1.00f);  // Dark background for headers (matches MenuBarBg/TitleBg)
    colors[ImGuiCol_TableBorderStrong] = ImVec4(0.30f, 0.30f, 0.35f, 1.00f);  // Outer borders (thicker)
    colors[ImGuiCol_TableBorderLight]  = ImVec4(0.22f, 0.22f, 0.27f, 1.00f);  // Inner borders (thinner)
    colors[ImGuiCol_TableRowBg]        = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);  // Transparent (use WindowBg)
    colors[ImGuiCol_TableRowBgAlt]     = ImVec4(1.00f, 1.00f, 1.00f, 0.02f);  // Very subtle white tint for alternating rows


    // Style tweaks
    style.WindowRounding = 0.0f;
    // style.WindowRounding = 5.0f;
    style.FrameRounding = 0.0f;
    // style.FrameRounding = 5.0f;
    style.GrabRounding = 5.0f;
    style.TabRounding = 0.0f;
    style.PopupRounding = 5.0f;
    style.ScrollbarRounding = 5.0f;
    style.WindowPadding = ImVec2(10, 10);
    style.FramePadding = ImVec2(6, 4);
    style.ItemSpacing = ImVec2(8, 6);
    style.WindowBorderSize = 1.0f;
    style.PopupBorderSize = 1.0f;
    style.ChildBorderSize = 1.0f;
    style.SeparatorTextBorderSize = 1.0f;  // Thinner separator line for SeparatorText
  }

  inline void SetLightTheme() {
    ImGuiStyle& style = ImGui::GetStyle();
    ImVec4* colors = style.Colors;

    ImVec4 bg = ImVec4(0.96f, 0.96f, 0.96f, 1.00f); // very light gray
    colors[ImGuiCol_WindowBg] = bg;

    style.WindowRounding = 5.0f;
    style.FrameRounding = 5.0f;
    style.GrabRounding = 5.0f;
    style.TabRounding = 0.0f;
    style.PopupRounding = 5.0f;
    style.ScrollbarRounding = 5.0f;
    style.WindowPadding = ImVec2(10, 10);
    style.FramePadding = ImVec2(6, 4);
    style.ItemSpacing = ImVec2(8, 6);
    style.WindowBorderSize = 1.0f;
    style.PopupBorderSize = 1.0f;
    style.ChildBorderSize = 1.0f;
    style.SeparatorTextBorderSize = 1.0f;  // Thinner separator line for SeparatorText

    ImVec4 bg_dark = ImVec4(0.92f, 0.92f, 0.92f, 1.00f); // light gray
    ImVec4 bg_light = ImVec4(1.00f, 1.00f, 1.00f, 1.00f); // white
    ImVec4 text = ImVec4(0.15f, 0.15f, 0.15f, 1.00f); // almost black
    ImVec4 text_muted = ImVec4(0.4f, 0.4f, 0.4f, 1.00f);    // mid-gray
    ImVec4 border = ImVec4(0.6f, 0.6f, 0.6f, 1.00f);    // light gray
    ImVec4 border_muted = ImVec4(0.7f, 0.7f, 0.7f, 1.00f);    // lighter gray
    ImVec4 primary = ImVec4(0.55f, 0.70f, 0.95f, 1.00f); // lighter office blue
    ImVec4 secondary = ImVec4(0.54f, 0.41f, 0.41f, 1.00f); // oklch(0.4 0.1 62)
    ImVec4 danger = ImVec4(0.50f, 0.35f, 0.30f, 1.00f); // oklch(0.5 0.05 30)
    ImVec4 warning = ImVec4(0.50f, 0.50f, 0.30f, 1.00f); // oklch(0.5 0.05 100)
    ImVec4 success = ImVec4(0.30f, 0.50f, 0.30f, 1.00f); // oklch(0.5 0.05 160)
    ImVec4 info = ImVec4(0.40f, 0.40f, 0.50f, 1.00f); // oklch(0.5 0.05 260)

    // Apply to ImGui color slots
    colors[ImGuiCol_Text] = text;
    colors[ImGuiCol_TextDisabled] = text_muted;
    colors[ImGuiCol_InputTextCursor] = text; // Match text color
    colors[ImGuiCol_WindowBg] = bg;
    colors[ImGuiCol_ChildBg] = bg;
    colors[ImGuiCol_PopupBg] = bg_light;
    colors[ImGuiCol_Border] = border;
    colors[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
    colors[ImGuiCol_FrameBg] = bg_dark;
    colors[ImGuiCol_FrameBgHovered] = border_muted;
    colors[ImGuiCol_FrameBgActive] = border;
    colors[ImGuiCol_TitleBg] = bg_dark;
    colors[ImGuiCol_TitleBgActive] = bg_dark;
    colors[ImGuiCol_TitleBgCollapsed] = bg_dark;
    colors[ImGuiCol_MenuBarBg] = bg_dark; // Slightly darker than window bg for contrast
    colors[ImGuiCol_ScrollbarBg] = bg_dark;
    colors[ImGuiCol_ScrollbarGrab] = border;
    colors[ImGuiCol_ScrollbarGrabHovered] = border_muted;
    colors[ImGuiCol_ScrollbarGrabActive] = text_muted;
    colors[ImGuiCol_CheckMark] = primary;
    colors[ImGuiCol_SliderGrab] = primary;
    colors[ImGuiCol_PlotHistogram] = primary;
    colors[ImGuiCol_SliderGrabActive] = primary;
    colors[ImGuiCol_Button] = bg_dark;
    colors[ImGuiCol_ButtonHovered] = border_muted;
    colors[ImGuiCol_ButtonActive] = border;
    colors[ImGuiCol_Header] = primary;
    colors[ImGuiCol_HeaderHovered] = primary;
    colors[ImGuiCol_HeaderActive] = primary;
    colors[ImGuiCol_Separator] = border;
    colors[ImGuiCol_SeparatorHovered] = border_muted;
    colors[ImGuiCol_SeparatorActive] = border;
    colors[ImGuiCol_ResizeGrip] = border;
    colors[ImGuiCol_ResizeGripHovered] = border_muted;
    colors[ImGuiCol_ResizeGripActive] = text;
    colors[ImGuiCol_Tab] = bg_dark;
    colors[ImGuiCol_TabHovered] = primary;
    colors[ImGuiCol_TabActive] = primary;
    colors[ImGuiCol_TabUnfocused] = bg_dark;
    colors[ImGuiCol_TabUnfocusedActive] = border;
    colors[ImGuiCol_DockingPreview] = primary;
    colors[ImGuiCol_DockingEmptyBg] = bg;
    colors[ImGuiCol_PlotLines] = secondary;
    colors[ImGuiCol_PlotLinesHovered] = danger;
    colors[ImGuiCol_TextSelectedBg] = primary;
    colors[ImGuiCol_DragDropTarget] = info;
    colors[ImGuiCol_NavHighlight] = primary;
    colors[ImGuiCol_NavWindowingHighlight] = ImVec4(1, 1, 1, 0.7f);
    colors[ImGuiCol_NavWindowingDimBg] = ImVec4(0.8f, 0.8f, 0.8f, 0.2f);
    colors[ImGuiCol_ModalWindowDimBg] = ImVec4(0.8f, 0.8f, 0.8f, 0.35f);
    colors[ImGuiCol_TableHeaderBg]     = bg_dark;
    colors[ImGuiCol_TableBorderStrong] = border; 
    colors[ImGuiCol_TableBorderLight]  = border_muted;
    style.SeparatorTextBorderSize = 1.0f;  // Thinner separator line for SeparatorText
  }

} // end Themes namespace
